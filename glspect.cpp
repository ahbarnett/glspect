// glSpect - OpenGL real-time spectrogram
// Alex Barnett 12/21/10, based on Luke Campagnola's nice 2005 glScope.
// ahb@math.dartmouth.edu
// Tweaked 10/13/11 for inverse video, etc.
// 1/25/14 fix -EPIPE snd_pcm_readi() error. nchannels=2, winf=0 case, etc
//        added color_byte to clean up color mapping; color map.
// 1/17/16: freq indicator line via left-button

/* Notes: be sure to set vSync wait in graphics card (eg NVIDIA) OpenGL settings
 */

/* ISSUES: (2011)
 * occasional dropped audio every few secs - why? (1470 vs 1472? issue)
 * jitter in signal graph scrolling
 * GlutGameMode isnt' setting refresh to 60Hz, rather 75Hz.
 * Better than glDrawPixels (which sends all data every frame to GPU) woudl be
    to pass the data as a texture and scroll it in the GPU (a convolution?),
    modifying only one row each time. This would be super low CPU usage!
 * Use GL_ARB_pixel_buffer_object for fast glDrawPixels or textures via DMA?
 * Add color?
 * add playback of audio file, or jack into audio playback?
 * glDrawpixels is deprecated in OpenGl >3.0. THat's annoying.
   Eg: https://www.opengl.org/discussion_boards/showthread.php/181907-drawing-individual-pixels-with-opengl
  "The modern way to do it is to store the data in a texture then draw a pair of textured triangles (quads are also deprecated)."
 */

/* SOLVED ISSUES: (2011)
 * Creation and initialization of global scn happens before main(), bad, since
     couldn't set t_memory in cmd line! Ans: use trivial creator/destructor,
     and call other init/close routines from main().
 */

#include <GL/glut.h>
#include <GL/gl.h>
#include <alsa/asoundlib.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include <sys/time.h>
#include <math.h>
#include <fftw3.h>

#define PI              3.14159265358979323846

class audioInput;      // declarations (things needed before defined)
void computespecslice(audioInput *ai);
int chooseTics(float lo, float range, float fac, float *tics);

static int verb;         // global verbosity
struct Param {           // parameters object
  int windowtype;        // DFT windowing function type
  int twowinsize;        // power of 2 giving DFT win_size (N) in samples
};                       // note trailing ;
Param param;           // global parameters object

const char *notenames[] = { "C","C#","D","Eb","E","F","F#","G","G#","A","Bb","B"};  // 1/17/16

float timeDiff(timeval a, timeval b) {
  return (float)(b.tv_sec - a.tv_sec) + (float)(b.tv_usec - a.tv_usec) * 0.000001F;
}

void drawText(float x, float y, char *string) {
  int len, i;
  glRasterPos2f(x, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++)
  { 
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
  }
}

void smallText(float x, float y, char *string) {
  int len, i;
  glRasterPos2f(x, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++)
  {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, string[i]);
  }
}

///////////// AUDIOINPUT class (no references to scn) ////////////////////
class audioInput {
  public:
  char *chunk;
  float *b, *specslice, *bwin, *winf, *sg;
  char *sgb;
  int b_ind, b_size, n_f, n_tw, sg_size, win_size;
  float dt, t_memory, Hz_per_pixel;
  
  bool quit, pause;
  pthread_t capture_thread;
  fftwf_plan fftw_p;

  // Audio device data (modified from ALSA tutorial)
  int bytes_per_frame, frames_per_period, nperiods, channels;
  int req_rate, rate;   /* Requested and actual sample rate */
  int dir;          /* rate == req_rate --> dir = 0 */
                    /* rate < req_rate  --> dir = -1 */
                    /* rate > req_rate  --> dir = 1 */
  snd_pcm_uframes_t period_size;     // Period size (bytes)
  snd_pcm_uframes_t req_size, size;  // requested and actual ALSA buffer size
  snd_pcm_t *pcm_handle;        /* Handle for the PCM device */ 
  snd_pcm_stream_t stream;     /* Playback stream */
  /* This structure contains information about    */
  /* the hardware and can be used to specify the  */      
  /* configuration to be used for the PCM stream. */ 
  snd_pcm_hw_params_t *hwparams;
  /* Name of the PCM device, like plughw:0,0          */
  /* The first number is the number of the soundcard, */
  /* the second number is the number of the device.   */
  char *pcm_name;

  void setupWindowFunc(float *w, int N) {
    float W;
    int i;
    if (verb) printf("windowtype=%d\n", param.windowtype);
    switch (param.windowtype) {
    case 0:  // no window (crappy frequency spillover)
      for( i=0; i<N; ++i)
	w[i] = 1.0F;
      break;
    case 1:  // Hann window (C^1 cont, so third-order tails)
      W = N/2.0F;
      for( i=0; i<N; ++i)
	w[i] = (1.0F + cos(PI*(i-W)/W))/2;
      break;
    case 2:  // truncated Gaussian window (Gaussian tails + exp small error)
      W = N/5.0F;    // width: keep small truncation but wide to not waste FFT
      for( i=0; i<N; ++i) w[i] = exp(-(i-N/2)*(i-N/2)/(2*W*W));
      break;
    default:
      fprintf(stderr, "unknown windowtype!\n");
    }
  }

  audioInput() {            // nontrivial constructor
    quit = false;
    pause = false;
    channels = 2;       // Had to change to stereo for System76 ! (was mono)
    bytes_per_frame = 2 * channels;      // 16-bit
    req_rate = 44100;         // audio sampling rate in Hz
    frames_per_period = (int)(req_rate/60.0);   // 735 = 44100Hz/60fps assumed
    nperiods = 2;             // >=2, see ALSA manual
    t_memory = 20.0;            // memory of our circular buffer in secs

    period_size = frames_per_period * bytes_per_frame;
    chunk = new char[period_size];   // raw data buffer for PCM read: 1 period
    req_size = frames_per_period * nperiods; // ALSA device buffer size (frames)
    b_ind = 0;                // integer index where to write to in buffer    
    if( initDevice() < 0 ) // set up sound card for recording (sets rate, size)
      exit(1);
    dt = 1.0 / (float)rate;   // sampling period
    b_size = (int)(t_memory * rate);   // buffer size
    if (verb) printf("memory buffer size %d samples\n", b_size);
    b = new float[b_size];    // our circular audio buffer
    
    win_size = 1<<param.twowinsize;     // FFT size
    bwin = new float[win_size];         // windowed recent audio data
    winf = new float[win_size];     
    setupWindowFunc(winf, win_size);    // windowing function
     // set up fast in-place single-precision real-to-half-complex DFT:
    fftw_p = fftwf_plan_r2r_1d(win_size, bwin, bwin, FFTW_R2HC, FFTW_MEASURE);
    n_f = 560;  // # freqs   ...spectrogram stuff
    specslice = new float[n_f];
    n_tw = 940; // # time windows: should be multiple of 4 for glDrawPixels
    sg_size = n_f * n_tw;
    sg = new float[sg_size];  // spectrogram array float
    sgb = new char[sg_size];  // spectrogram array 8-bit
    //for( int i=0; i<sg_size; ++i )  // fill with random for now
    //sg[i] = (float)(random()/(float)RAND_MAX);
    //for( int i=0; i<sg_size; ++i )  // fill with random for now
    //sgb[i] = (int)(256.0*random()/(float)RAND_MAX);
    Hz_per_pixel = 1.0F / (win_size*dt);
    if (verb) printf("Hz per pixel = %.3f\n", Hz_per_pixel);
 
    // Start recording thread... runs independently, writing data into ai->b
    pthread_create(&capture_thread, NULL, audioCapture, (void*)this); // this?
  }
  ~audioInput() {             // destructor
    snd_pcm_close (pcm_handle);
    fftwf_destroy_plan(fftw_p);
  }
  
  int initDevice() {  // ........ set up sound card for recording ........
    // ALSA tutorial, taken from http://www.suse.de/~mana/alsa090_howto.html
    
    stream = SND_PCM_STREAM_CAPTURE;
    /* Init pcm_name. Of course, later you */
    /* will make this configurable ;-)     */
    pcm_name = strdup("plughw:0,0");
    /* Allocate the snd_pcm_hw_params_t structure on the stack. */
    snd_pcm_hw_params_alloca(&hwparams);
    /* Open PCM. The last parameter of this function is the mode. */
    /* If this is set to 0, the standard mode is used. Possible   */
    /* other values are SND_PCM_NONBLOCK and SND_PCM_ASYNC.       */ 
    /* If SND_PCM_NONBLOCK is used, read / write access to the    */
    /* PCM device will return immediately. If SND_PCM_ASYNC is    */
    /* specified, SIGIO will be emitted whenever a period has     */
    /* been completely processed by the soundcard.                */
    if (snd_pcm_open(&pcm_handle, pcm_name, stream, 0) < 0) {
      fprintf(stderr, "Error opening PCM device %s\n", pcm_name);
      return(-1);
    }
    /* Init hwparams with full configuration space */
    if (snd_pcm_hw_params_any(pcm_handle, hwparams) < 0) {
      fprintf(stderr, "Can not configure this PCM device.\n");
      return(-1);
    }
    /* Set access type. This can be either    */
    /* SND_PCM_ACCESS_RW_INTERLEAVED or       */
    /* SND_PCM_ACCESS_RW_NONINTERLEAVED.      */
    /* There are also access types for MMAPed */
    /* access, but this is beyond the scope   */
    /* of this introduction.                  */
    if (snd_pcm_hw_params_set_access(pcm_handle, hwparams, SND_PCM_ACCESS_RW_INTERLEAVED) < 0) {
      fprintf(stderr, "Error setting access.\n");
      return(-1);
    }
    
    /* Set sample format */
    if (snd_pcm_hw_params_set_format(pcm_handle, hwparams, SND_PCM_FORMAT_S16_LE) < 0) {
      fprintf(stderr, "Error setting format.\n");
      return(-1);
    }
    
    /* Set sample rate. If the requested rate is not supported */
    /* by the hardware, use nearest possible rate.         */ 
    rate = req_rate;
    if (snd_pcm_hw_params_set_rate_near(pcm_handle, hwparams, (uint*)&rate, 0) < 0) {
      fprintf(stderr, "Error setting rate.\n");
      return(-1);
    }
    if (rate != req_rate) {
      fprintf(stderr, "The rate %d Hz is not supported by your hardware.\n \
                        ==> Using %d Hz instead.\n", req_rate, rate);
    }
    
    /* Set number of channels */
    if (snd_pcm_hw_params_set_channels(pcm_handle, hwparams, channels) < 0) {
      fprintf(stderr, "Error setting channels.\n");
      return(-1);
    }
    
    /* Set number of periods. Periods used to be called fragments. */ 
    if (snd_pcm_hw_params_set_periods(pcm_handle, hwparams, nperiods, 0) < 0) {
      fprintf(stderr, "Error setting number of periods.\n");
      return(-1);
    }
    /* Set buffer size (in frames). The resulting latency is given by */
    /* latency = period_size * nperiods / (rate * bytes_per_frame)     */
    size = req_size;
    if (snd_pcm_hw_params_set_buffer_size_near(pcm_handle, hwparams, &size) < 0) {
      fprintf(stderr, "Error setting buffersize.\n");
      return(-1);
    }
    if( size != req_size ) {
      fprintf(stderr, "Buffer size %d is not supported, using %d instead.\n", (int)req_size, (int)size);
    }
    
    /* Apply HW parameter settings to PCM device and prepare device  */
    if (snd_pcm_hw_params(pcm_handle, hwparams) < 0) {
      fprintf(stderr, "Error setting HW params.\n");
      return(-1);
    }
    return 1;
  } // ........................................
  
  void quitNow() {
    quit = true;
    //      pthread_kill_other_threads_np();
    snd_pcm_close (pcm_handle);
  }
  
  union byte {                    // used to convert from signed to unsigned
    unsigned char uchar_val;
    char char_val;
  };
  
  int mod( int i ) {  // true modulo (handles negative) into our buffer b
    // wraps i to lie in [0, (b_size-1)]. rewritten Barnett
    int r = i % b_size;
    if (r<0)
      r += b_size;
    return r;
  }
  
  static void* audioCapture(void* a) { //-------- capture: thread runs indep --
    // still mostly Luke's code, some names changed. Aims to read 1 "period"
    // (ALSA device setting) into the current write index of our ai->b buffer.
    fprintf(stderr, "audioCapture thread started...\n");
    audioInput* ai = (audioInput*) a;  // shares data with main thread = cool!
    
    float inv256 = 1.0 / 256.0;
    float inv256_2 = inv256*inv256;
    
    while( ! ai->quit ) {  // loops around until state of ai kills it
      int n;
      if( ! ai->pause ) {
	// keep trying to get exactly 1 "period" of raw data from sound card...
	while((n = snd_pcm_readi(ai->pcm_handle, ai->chunk, ai->frames_per_period)) < 0 ) {
	  //	  if (n == -EPIPE) fprintf(stderr, "Overrun occurred: %d\n", n); // broken pipe
	  fprintf(stderr, "Error occured while recording: %s\n", snd_strerror(n));

	  //n = snd_pcm_recover(ai->pcm_handle, n, 0); // ahb

	  //fprintf(stderr, "Error occured while recording: %s\n", snd_strerror(n));
	  snd_pcm_prepare(ai->pcm_handle);
	  //fprintf(stderr, "Dropped audio data (frames read n=%d)\n", n);
	}  // n samples were got
	if (verb>1) printf("snd_pcm_readi got n=%d frames\n", n);

	byte by;
	int write_ptr, read_ptr;
	for( int i = 0; i < n; i++ ) { // read chunk into our buffer ai->b ...
	  read_ptr = i * ai->bytes_per_frame;
	  write_ptr = ai->mod(ai->b_ind + i); // wraps around
	  by.char_val = ai->chunk[read_ptr];
	  // compute float in [-1/2,1/2) from 16-bit raw... (LSB unsigned char)
	  ai->b[write_ptr] = (float)ai->chunk[read_ptr+1]*inv256 + (float)by.uchar_val*inv256_2;
	}
	ai->b_ind = ai->mod(ai->b_ind+n);  // update index (in one go)

	computespecslice(ai); // compute spectral slice of recent buffer history
      }
      else {
	usleep(10000);  // wait 0.01 sec if paused (keeps thread CPU usage low)
      }
    }
    fprintf(stderr, "audioCapture thread exiting.\n");
  }                          // ----------------------- end capture thread ----
};


//////////////////////////////////////// SCENE CLASS /////////////////////////
class scene
{
  public:
  float viewport_size[2];  // gl window width, height
  float mouse[3];
  bool pause;
  audioInput* ai;
  float color_scale[2];    // sg image log color mapping params
  bool diagnose;           // medical stuff
  char diagnosis[1000];
  time_t fps_tic;             // FPS time reference
  int frameCount, fps;       // frames per sec and counter
  int scroll_fac, scroll_count, colormode;
  float run_time;
  timeval start_time;
  int freqind;     // toggles frequency readoff line(s), either 0,1,2

  scene() { }        // trivial constructor, to prevent global init fiasco
  ~scene() { }      // destructor
  
  void init() {       // meaningful constructor - could move some to main to
                      // parse cmd line
    pause = false;
    color_scale[0] = 100.0;     // 8-bit intensity offset
    color_scale[1] = 255/120.0;     // 8-bit intensity slope (per dB units)
    ai = new audioInput();
    fps_tic = time(NULL); frameCount = fps = 0;
    colormode = 0;
    scroll_count = 0;
    gettimeofday(&start_time, NULL);
    run_time = 0.0;
    freqind = 0;
    diagnose = 0;                         // initialize Medical stuff
    strcpy(diagnosis, "no diagnosis ..."); 
  }
  
  void graph() // simple graph of ai->float_data ............................
  {
    float tshow = 0.1;    // only show the most recent tshow secs of samples
    glDisable(GL_DEPTH_TEST); glDisable(GL_BLEND);
    glDisable(GL_LINE_SMOOTH); glLineWidth(2);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix(); // use modelview matrix to xform [-tshow,0]x[-1,1] somewhere
    glTranslatef(0.97,0.1,0);
    glScalef(4.0,1.0,1.0);  // x-scale for time-units, y-scale is 1
    char xlab[] = "t(s)", ylab[] = "";   // axis labels
    //drawAxes(-tshow, 0, -0.1, 0.1, 1, 1, 0, 0, xlab, ylab);   // not needed?
    glColor4f(0.4, 1.0, 0.6, 1);
    glBegin(GL_LINE_STRIP);
    float x = -tshow;
    int ilo = ai->b_ind-(int)(tshow*ai->rate);
    int ihi = ai->b_ind; // NB get now since capture thread may change it!
    for( int i=ilo; i<ihi; i++ ) {   // graph the most recent piece of buffer
      glVertex2f(x, ai->b[ai->mod(i)]);
      x += ai->dt;   // x time increment
    }
    glEnd();
    glPopMatrix();
  }           // ............................................................

  void slice() // simple graph of spectral slice with log power scale
  {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix(); // use modelview matrix to xform [-tshow,0]x[-1,1] somewhere
    glTranslatef(0.15,0.11,0);   // decide the location of 0,0 of the graph
    float max_Hz = ai->Hz_per_pixel*ai->n_f;
    glScalef(0.3/max_Hz,0.0015,1.0);  // x-scale for freq-units, y-scale for dB
    glLineWidth(1);
    glColor4f(0.5, 0.4, 0.2, 1);    // lines showing spectrogram color range
    glBegin(GL_LINES);
    float dB_min = -color_scale[0]/color_scale[1];
    float dB_max = (255.0-color_scale[0])/color_scale[1];
    glVertex2f(0,dB_min);glVertex2f(max_Hz,dB_min);
    glVertex2f(0,dB_max);glVertex2f(max_Hz,dB_max);
    glEnd();
    glColor4f(1.0, 0.8, 0.3, 1);   // the graph
    glBegin(GL_LINE_STRIP);
    for( int i=0; i<ai->n_f; i++ )
      glVertex2f(i*ai->Hz_per_pixel, 20*log10(ai->specslice[i]));
    glEnd();
    char xlab[] = "f(Hz)", ylab[] = "I(dB)";   // axis labels
    drawAxes(0-1e-7, max_Hz, -50, 50, 1.0, 1.0, 1, 0, xlab, ylab);
    glPopMatrix();
  }

  void drawAxes(float x0, float x1, float y0, float y1, float xf, float yf, bool box, bool grid, char* xlab, char* ylab)  // draws axes in [x0,x1]x[y0,y1]
  // box: 0 coord axes no box, 1 left-bottom axes no box, 2 box no axes
  // xf, yf are density factors, roughly given by screen size of object.
  // xlab, ylab are labels
  {
    int n_tics, i, len;
    char label[20];
    float tics[100];
    float x_pix_size = 2.0*(x1-x0)/(viewport_size[0]*xf); // fudge fac
    float y_pix_size = 2.0*(y1-y0)/(viewport_size[1]*yf);

    glColor4f(0.7, 1.0, 1.0, 1);
    glDisable(GL_LINE_SMOOTH); glLineWidth(1);
    if (box) { glBegin(GL_LINE_STRIP);  // box and LB axes variants
      glVertex2f(x0,y1); glVertex2f(x0,y0); glVertex2f(x1,y0);
      if (box>1) { glVertex2f(x1,y1); glVertex2f(x0,y1); }
    glEnd();
    } else { glBegin(GL_LINES); // coord axes lines
    glVertex2f(0,y0); glVertex2f(0,y1); glVertex2f(x0,0); glVertex2f(x1,0);
    glEnd(); }
    smallText(x1+8*x_pix_size,y0-8*y_pix_size,xlab); // axes labels
    smallText(x0-8*x_pix_size,y1+8*y_pix_size,ylab);
    
    n_tics = chooseTics(x0,x1-x0,xf,tics);  // x-axis tics
    float xtic_top = (box==0) ? 0.0 : y0;
    float xtic_bot = xtic_top - 0.02*(y1-y0)/yf;
    glBegin(GL_LINES);
    for (i=0;i<n_tics;++i) {
      glVertex2d(tics[i], xtic_bot); glVertex2d(tics[i], xtic_top);
    }
    glEnd();
    for (i=0;i<n_tics;++i) {
      sprintf(label,"%.6g",tics[i]);
      smallText((float)(tics[i] - 4*strlen(label)*x_pix_size), \
	     (float)(xtic_bot - 12*y_pix_size), label);
    }

    n_tics = chooseTics(y0,y1-y0,yf,tics);  // y-axis tics
    float ytic_top = (box==0) ? 0.0 : x0;
    float ytic_bot = ytic_top - 0.02*(x1-x0)/xf;
    glBegin(GL_LINES);
    for (i=0;i<n_tics;++i) {
      glVertex2d(ytic_bot, tics[i]); glVertex2d(ytic_top, tics[i]);
    }
    glEnd();
    for (i=0;i<n_tics;++i) {
      sprintf(label,"%.6g",tics[i]);
      smallText((float)(ytic_bot - 8*strlen(label)*x_pix_size), \
	     (float)(tics[i] - 4*y_pix_size), label);
    }
  }

  void spectrogram() // show spectrogram as 2D pixel array, w/ axes............
  {
    float x0=0.05, y0=0.22;    // bot-left location in viewport (as unit square)
    float currfreq, linefreq;  char str[50];  // for freqind
    int nummults, i, notenum, octave;   // for freqind
    glRasterPos2f(x0,y0);
    glDisable(GL_DEPTH_TEST); glDisable(GL_BLEND);
    //glPixelZoom(1.f,1.f); // zoom: non-integers are quite a bit slower (ie twice), and seems to have permanently worsened the DrawPixels speed by factor 2!
    //glPixelZoom(1.F,2.F);
    //float z; glGetFloatv(GL_ZOOM_Y, &z); printf("%f\n",z); // read a zoom
    if (colormode<2)   // b/w
      glDrawPixels(ai->n_tw, ai->n_f, GL_LUMINANCE, GL_UNSIGNED_BYTE, ai->sgb);
    else
      glDrawPixels(ai->n_tw, ai->n_f, GL_RGB, GL_UNSIGNED_BYTE_3_3_2, ai->sgb);

    //glDrawPixels(ai->n_tw, ai->n_f, GL_RED, GL_UNSIGNED_BYTE, ai->sgb); //test

    // from: http://www.gamedev.net/topic/510914-glcolortable-with-textures/
    //float c[256][3];
    //for (int i=0;i<256;i++) {
    //  c[i][0] = ((float)i)/256;
    //  c[i][1] = 1.0f - ((float)i)/256;
    //  c[i][2] = ((float)i)/256;
    //}
    //glColorTable(GL_COLOR_TABLE, GL_RGB, 256, GL_RGB, GL_FLOAT, c);
    //glEnable(GL_COLOR_TABLE);
    //glDrawPixels(ai->n_tw, ai->n_f, GL_BITMAP,GL_COLOR_INDEX, ai->sgb);
    // fails!

    glMatrixMode(GL_MODELVIEW); glPushMatrix(); // xform axes into right place
    float secs_per_pixel = scroll_fac / 60.0; // (float)fps, or longer FPS mean?
    float max_Hz = ai->Hz_per_pixel*ai->n_f;
    float max_t = secs_per_pixel*ai->n_tw;
    glTranslatef(x0 + (ai->n_tw-run_time/secs_per_pixel)/(float)viewport_size[0], y0, 0);
    glScalef(1.0F/(viewport_size[0]*secs_per_pixel),1.0F/(viewport_size[1]*ai->Hz_per_pixel),1); // makes units match the spectrogram t (scrolling) & freq axes
    // plots now occur in physical t (s) and f (Hz) units, relative to start_t..
    char xlab[] = "t(s)", ylab[] = "f(Hz)";
    drawAxes(run_time-max_t, run_time, 0 - 1e-7, max_Hz, 2.0, 2.0, 2, 0, xlab, ylab);
    // 1e-7 is hack to show 0 Hz
    if (freqind) {        // horiz freq readout line, in Hz/sec coords, 1/17/16
      // reverse-engineer the freq from current mouse y in pixels (yuk)...
      currfreq = ai->Hz_per_pixel*(viewport_size[1]*(1-y0)-mouse[1]);
      if (currfreq>0.0) {       // only show if meaningful freq
	nummults = (freqind>1) ? 10 : 1;
	for (i=1;i<=nummults;++i) {      // do either one or many lines
	  linefreq = i*currfreq;
	  if (i==1)
	    glColor4f(0.8, 0.5, 0.0, 1); // orange
	  else
	    glColor4f(0.5, 0.3, 0.0, 1); // darker orange
	  glDisable(GL_LINE_SMOOTH); glLineWidth(1);
	  glBegin(GL_LINES);
	  glVertex2f(run_time-max_t, linefreq);
	  glVertex2f(run_time, linefreq);
	  glEnd();
	  glEnable(GL_BLEND); // text label: overwrite entirely
	  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
	  notenum = (int)roundf(12 * logf(linefreq/261.626)/logf(2.0));
	  octave = 4 + notenum/12; // relies on int division
	  if (octave<0) octave = (int)NAN;
	  sprintf(str, "  %.d   %s%d",(int)roundf(linefreq),notenames[(notenum+1200) % 12],octave);  // use global lookup table of note names
	  smallText(run_time-max_t,linefreq+3.0*ai->Hz_per_pixel,str);
	}
      }
    }
    glPopMatrix();
  }               // .......................................................

  char color_byte(float x)  // convert spectrogram entry to 8bit color char
  // barnett 1/25/15
  {
    float fac = 20.0 * color_scale[1]; // color_scale[1] units: intensity/dB
    int k = (int)(color_scale[0] + fac*log10(x));
    if (k>255) k=255; else if (k<0) k=0;
    if (colormode==0)           // b/w
      return (char)k;
    else if (colormode==1)      // inverse b/w
      return (char)(255-k);
    else if (colormode==2) {    // color (pack into 3_3_2 RGB format) .. SLOW?
      float a = k/255.0f;       // 0<a<1. now map from [0,1] to rgb in [0,1]
      float r = 5*(a-0.2); if (r<0) r=0.0; else if (r>=1) r=.955; // clip
      float g = 5*(a-0.6); if (g<0) g=0.0; else if (g>=1) g=.995;
      float b = 5*a;
      if (a>0.8) b = 5*(a-0.8); else if (a>0.4) b = 5*(0.6-a);      
      if (b<0) b=0.0; else if (b>=1) b=.995;
      return (char)(b*4 + 4*((int)(g*8)) + 32*((int)(r*8)));  // pack to 8bits
    }
  }

  void regen_sgb()   // recompute 8bit spectrogram from float spectrogram
  // can only be done at a few FPS, so only do it when color scale changes
  {
    int i, j, n = ai->n_tw;
    for ( i=0; i<n; ++i)  // loop over whole sg array
      for ( j=0; j<ai->n_f; ++j)
	ai->sgb[j*n + i] = color_byte(ai->sg[j*n + i]);
  }

}; //////////////////////////////////////////// END SCENE CLASS ////////////////

static scene scn;    // global scene which contains everything (eg via scn.ai)

int chooseTics(float lo, float range, float fac, float *tics)
/* returns the number of tics, and their locations in zero-indexed tics[].
 * Barnett 99/10/20, single-prec version of ~visu/viewer/viewer.c
 * But, has density fudge factor fac (if 0, chooses default value).
 */
{
    int i, n_tics, tic_start;
    float exponent, logr, spacing;
    if (fac==0.0) fac = 1.0;

    /* adjust the range multiplier here to give good tic density... */
    logr = log10(range * 0.4 / fac);
    exponent = floor(logr);
    spacing = pow(10.0, exponent);
    if (logr-exponent > log10(5.0))
        spacing *= 5.0;
    else if (logr-exponent > log10(2.0))
        spacing *= 2.0;
    
    /* (int) and copysign trick is to convert the floor val to an int... */
    tic_start = (int)(copysign(0.5,lo) + 1.0 + floor(lo / spacing));
    
    n_tics = (int)(1.0 + (lo + range - tic_start*spacing)/spacing);
    for (i=0;i<n_tics;++i) {
        tics[i] = spacing * (tic_start + i);
    }
    return n_tics;
}


void drawText()       // Medical text and other text overlays................
{
  glEnable (GL_BLEND); // glBlendFunc (GL_ONE,GL_ZERO); // overwrite entirely
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

  glColor4f(1, .4, .4, 1.0); // write FPS, other params... coords in unit square
  char str[50]; sprintf(str, "%d FPS",scn.fps); smallText(0.92,0.96,str);
  sprintf(str, "gain offset %.0f dB",scn.color_scale[0]);
  smallText(0.02,0.04,str);
  sprintf(str, "dyn range  %.1f dB",255.0/scn.color_scale[1]);
  smallText(0.02,0.02,str);

  if (scn.diagnose) {  // show diagnosis
    glColor4f(.2, .2, .2, 0.8);  // transparent box
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix(); // use modelview matrix to xform unit square to text box
      glTranslatef(0.5,0.5,0); glScalef(0.5,0.4,1); glTranslatef(-0.5,-0.5,0); 
      glBegin(GL_POLYGON);
      glVertex2f(0,0); glVertex2f(1,0); glVertex2f(1,1); glVertex2f(0,1); glVertex2f(0,0); // unit square
      glEnd();
      glColor4f(1, 1, 1, 1);                  // text
      drawText(0.2, 0.5, scn.diagnosis); // coords relative to box as unit sq
    glPopMatrix();
  }
}
 
void display()  //  GLUT's display routine: 2D layering by write order
{
  glClear(GL_COLOR_BUFFER_BIT); // no depth buffer
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, 1, 0, 1, -1, 1);  // l r b t n f
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  scn.graph();      // simple graph of signal
  scn.slice();   // spectral slice
  scn.spectrogram();   // note, overlays the signal and spectrum graphs
  drawText();    // overlaid
  // glFinish();   // wait for all gl commands to complete

  glutSwapBuffers(); // for this to WAIT for vSync, need enable in NVIDIA OpenGL
}

void computespecslice(audioInput *ai)   // windowed spectral slice!
  {
    int N = ai->win_size;             // transform length
    int nf = ai->n_f;              // # freqs to fill in powerspec
   
    if (0)
      for ( int i=0; i<nf; ++i )     // dummy just copy last n_f samples
	ai->specslice[i] = 100.0 * ai->b[ai->mod(ai->b_ind - nf + i)];

    else {                      // actual windowed power spectrum
      for ( int i=0; i<N; ++i )     // copy last N samples & mult by window func
	ai->bwin[i] = ai->winf[i] * ai->b[ai->mod(ai->b_ind - N + i)];

      fftwf_execute(ai->fftw_p); // do the already-planned fft

      if (nf>N/2) {fprintf(stderr,"window too short cf n_f!\n"); return;}
      ai->specslice[0] = ai->bwin[0]*ai->bwin[0];  // zero-freq has no imag
      for ( int i=1; i<nf; ++i )     // compute power spectrum from hc dft
	ai->specslice[i] = ai->bwin[i]*ai->bwin[i] + ai->bwin[N-i]*ai->bwin[N-i];
    }
  }

void scroll_sg(audioInput *ai, float *newcol)
// Scroll spectrogram data 1 pixel in t direction, add new col & make 8bit
{
  int i, j, n = ai->n_tw;

  for ( i=0; i<n-1; ++i)    // float: NB x (ie t) is the fast storage direc
    for ( j=0; j<ai->n_f; ++j)
      ai->sg[j*n + i] = ai->sg[j*n + i + 1];

  for ( i=0; i<n; ++i)    // scroll the 8bit char data too
    for ( j=0; j<ai->n_f; ++j)
      ai->sgb[j*n + i] = ai->sgb[j*n + i + 1];
  
  for ( j=0; j<ai->n_f; ++j)   // set the last col of float data
    ai->sg[j*n + n - 1] = newcol[j];
  
  for ( j=0; j<ai->n_f; ++j) {        // color xform for last col of 8bit data
    ai->sgb[j*n + n - 1] = scn.color_byte(newcol[j]);
  }
}

// --------------------------- IDLE FUNC ------------------------------------
void idle(void)   //    put numerical intensive part here? only scrolling
{
  ++scn.frameCount;        // compute FPS...
  time_t now = time(NULL);     // integer in seconds
  if (now>scn.fps_tic) {    // advanced by at least 1 sec?
    scn.fps = scn.frameCount;
    scn.frameCount = 0;
    scn.fps_tic = now;
  }

  if (!scn.ai->pause) {  // scroll every scroll_fac vSyncs, if not paused:
    timeval nowe;  // update runtime
    gettimeofday(&nowe, NULL);
    //scn.run_time = timeDiff(scn.start_time, nowe);  // update run_time data
    scn.run_time += 1/60.0;        // is less jittery, approx true
    scn.run_time = fmod(scn.run_time, 100.0);    // wrap time label after 100 secs

    ++scn.scroll_count;
    if (scn.scroll_count==scn.scroll_fac) {
      scroll_sg(scn.ai, scn.ai->specslice);   // add spec slice to sg & scroll
      scn.scroll_count = 0;
    }
  }

  glutPostRedisplay();  // trigger GLUT display func
}


void keyboard(unsigned char key, int xPos, int yPos)
{
  if( key == 27 || key == 'q') {  // esc or q to quit
    scn.ai->quitNow();
    exit(0);
  } else if( key == ' ' ) {
    scn.pause = ! scn.pause;
    scn.ai->pause = ! scn.ai->pause;
  } else if( key == 'd' ) {
    scn.diagnose = ! scn.diagnose;         // toggle text overlay
  } else if( key == ']' ) {               // speed up scroll rate
    if (scn.scroll_fac>1) {scn.scroll_fac--; scn.scroll_count=0; }
  } else if( key == '[' ) {
    if (scn.scroll_fac<50) scn.scroll_fac++;
  } else if( key == 'i' ) {
    scn.colormode = (scn.colormode + 1) % 3;     // spectrogram color scheme
    scn.regen_sgb();
  } else {
    fprintf(stderr, "pressed key %d\n", (int)key);
  }
}

void special(int key, int xPos, int yPos) {  // special key handling
  if( key == 102 ) { // rt
    scn.color_scale[1] *= 1.5; scn.regen_sgb(); // contrast
  } else if( key == 100 ) { // lt
    scn.color_scale[1] /= 1.5; scn.regen_sgb(); // contrast
  } else if( key == 103 ) { // dn
    scn.color_scale[0] -= 20; scn.regen_sgb();  // brightness
  } else if( key == 101 ) { // up
    scn.color_scale[0] += 20; scn.regen_sgb();  // brightness
  } else {
    fprintf(stderr, "pressed special key %d\n", key);
  }
}

void reshape(int w, int h) // obsolete in game mode
{
  scn.viewport_size[0] = w;
  scn.viewport_size[1] = h;
  if (verb) fprintf(stderr, "Setting w=%d, h=%d\n", w, h);
  glViewport(0, 0, w, h);
}

void mouse(int button, int state, int x, int y)
{
  scn.mouse[0] = x;  // this seems to be needed for correct panning
  scn.mouse[1] = y;
  scn.mouse[2] = button;
  if (verb) printf("Mouse event: button: %d  state: %d  pos: %d, %d\n", button, state, x, y);
  if (button==GLUT_LEFT_BUTTON && state==GLUT_DOWN) scn.freqind = 1; // toggle
  if (button==GLUT_LEFT_BUTTON && state==GLUT_UP) scn.freqind = 0;
  if (button==GLUT_RIGHT_BUTTON && state==GLUT_DOWN) scn.freqind = 2; // toggle w/ freq multiples shown
  if (button==GLUT_RIGHT_BUTTON && state==GLUT_UP) scn.freqind = 0;
  if (button==GLUT_MIDDLE_BUTTON && state==GLUT_UP) scn.regen_sgb();
}

void motion(int x, int y) // handles mouse drag effects
{
  int dx = (int)(x-scn.mouse[0]);
  int dy = (int)(y-scn.mouse[1]);
  if( scn.mouse[2] == GLUT_MIDDLE_BUTTON ) {   // controls color scale
    scn.color_scale[0] += dx/5.0;  // brightness
    scn.color_scale[1] *= exp(-dy/200.0); // contrast
  }
  scn.mouse[0] = x;
  scn.mouse[1] = y;
}

const char *helptext[] = {
  " glSpect: real-time OpenGL spectrogram.  Alex Barnett, Dec 2010\n",
  "                                         (based on Luke Campagnola glScope)\n\n",
  "Usage: glspect  [-f] [-v] [-sf <scroll_factor>] [-w <windowtype>] [-t twowinsize]\n\n",
  "Command line arguments:\n",
  "windowtype = \t0 (no window)\n\t\t1 (Hann)\n\t\t2 (Gaussian trunc at +-4sigma) (default)\n",
  "scroll_factor = 1,2,... # vSyncs (60Hz) to wait per scroll pixel (default 1)\n",
  "twowinsize = 11,12,...,16 is power of 2 giving FFT win_size N (default 13)\n\t(Note: this controls the vertical frequency resolution and range)\n\n", 
  "Keys & mouse: \tarrows or middle button drag - brightness/contrast\n",
  "\t\tleft button shows horizontal frequency readoff line\n",
  "\t\tright button shows horizontal frequency readoff with multiples\n",
  "\t\ti - cycles through color maps (B/W, inverse B/W, color)\n",
  "\t\tq or Esc - quit\n",
  "\t\t[ and ] - control horizontal scroll factor (rate)\n",
  NULL };

// ===========================================================================
int main(int argc, char** argv)
{
  int scrnmode = 0;  // 0 for window, 1 fullscreen       **Defaults go here**
  verb = 0;          // 0 silent, 1 debug, etc
  scn.scroll_fac = 2;    // how many vSyncs to wait before scrolling sg
  param.windowtype = 2;  // Gaussian
  param.twowinsize = 13; // 8192 samples (around 0.19 sec). Remains fixed

  for (int i=1; i<argc; ++i) {  // .....Parse cmd line options....
    if (!strcmp(argv[i], "-f"))  // option -f makes full screen
      scrnmode = 1;
    else if (!strcmp(argv[i], "-v"))  // option -v makes verbose
      verb = 1;
    else if (!strcmp(argv[i], "-sf")) {
      sscanf(argv[++i], "%d", &scn.scroll_fac);  // read in scroll factor
      if (scn.scroll_fac<1) scn.scroll_fac=1; // sanitize it
    } else if (!strcmp(argv[i], "-t")) {
      sscanf(argv[++i], "%d", &param.twowinsize);
      if (param.twowinsize<10) param.twowinsize=10; // sanitize it
      if (param.twowinsize>18) param.twowinsize=18;
    } else if (!strcmp(argv[i], "-w")) {
      sscanf(argv[++i], "%d", &param.windowtype);  // read in windowtype
      if (param.windowtype<0) param.windowtype=0; // sanitize it
    } else {            // misuse (or -h): print out usage text...
      fprintf(stderr, "bad command line option %s\n\n", argv[i]);
      for (int i=0; helptext[i]; i++) fprintf(stderr, "%s", helptext[i]);
      exit(1);
    }
  }
  scn.init();       // true constructor for global scn object

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  if (scrnmode==1) {  // implement desired screen mode
    //glutGameModeString( "1680x1050:24@60" ); // res, pixel depth, refresh
    glutGameModeString( "1024x768:24@60" ); // for XGA projector
    glutEnterGameMode();                     // start fullscreen game mode
  } else {
    glutInitWindowSize(1024, 768);  // window same size as XGA
    int mainWindow = glutCreateWindow("glSpect by Alex Barnett, Dec 2010");
    //  glutFullScreen();    // maximizes window, but is not game mode
  }

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(special);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutIdleFunc(idle);

  glutMainLoop();        // pass control to GLUT
  return 0;
}

// glscope, by Luke Canpagnola 2005, tweaked by Alex Barnett 2007-2015

#include <GL/glut.h>
#include <GL/gl.h>
#include <alsa/asoundlib.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include <sys/time.h>

 /*******************
   Helper functions
 ********************/

int sign(float val) {
  if( val < 0 ) 
    return -1;
  else
    return 1;
}

int rot(int a, int mod) {
  return (a+1) % mod;
}

float snap(float val, float grid) {
  float chop_val = (int)(val * 1000.0) / 1000.0; // workaround for small fp-error problems
  float snapped = (int)((chop_val/grid)+sign(chop_val) * 0.5) * grid;
  return(snapped);
}

float afloor_snap(float val, float grid) {
  float chop_val = (int)(val * 1000.0) / 1000.0; // workaround for small fp-error problems
  float snapped = (int)((chop_val/grid)) * grid;
  return(snapped);
}

float min(float a, float b) {
  return a>b ? b : a;
}
float max(float a, float b) {
  return a>b ? a : b;
}
typedef float Vector [4];

void doVertex(Vector v) {
  glVertex3f(v[0], v[1], v[2]);
}

float timeDiff(timeval a, timeval b) {
  return (float)(b.tv_sec - a.tv_sec) + ((float)(b.tv_usec - a.tv_usec) * .000001);
}

class audioInput {
  public:
    char* data;
    float* float_data;
    int buffer_blocks;
    int data_start;
    int data_size;
    int data_end;
    int data_write;
    int write_padding;
    int front_padding;
    
    bool quit;
    bool pause;
    pthread_t capture_thread;
    timeval last_write;
    float avg_write_interval;

      // Audio device data

    int bytes_per_frame;
    int read_frames, read_bytes;
    int rate; /* Sample rate */
    float dt;
    int exact_rate;   /* Sample rate returned by */
                      /* snd_pcm_hw_params_set_rate_near */ 
    int dir;          /* exact_rate == rate --> dir = 0 */
                      /* exact_rate < rate  --> dir = -1 */
                      /* exact_rate > rate  --> dir = 1 */
    int periods;       /* Number of periods */
    snd_pcm_uframes_t periodsize; /* Periodsize (bytes) */
    snd_pcm_uframes_t size, exact_size;  // buffer size   AHB changed from int
      /* Handle for the PCM device */ 
    snd_pcm_t *pcm_handle;
    /* Playback stream */
    snd_pcm_stream_t stream;
    /* This structure contains information about    */
    /* the hardware and can be used to specify the  */      
    /* configuration to be used for the PCM stream. */ 
    snd_pcm_hw_params_t *hwparams;
    /* Name of the PCM device, like plughw:0,0          */
    /* The first number is the number of the soundcard, */
    /* the second number is the number of the device.   */
    char *pcm_name;

    audioInput() {
      gettimeofday(&last_write, NULL);
      quit = false;
      pause = false;
      buffer_blocks = 50;
      bytes_per_frame = 2;
      read_frames = 1000;
      write_padding = read_frames * 4;
      front_padding = read_frames;
      avg_write_interval = 0.1;
      
      read_bytes = read_frames * bytes_per_frame;
      data = new char[read_bytes];
      
      data_start = read_frames * 2;
      data_end = data_size-1;
      data_write = 0;
      data_size = buffer_blocks * read_frames;
      float_data = new float[data_size];
      
      rate = 44100;
      dt = 1.0 / (float)rate;
      periods = 2;
      periodsize = 8192;
      size = exact_size = (periodsize * periods)>>2;
      
      if( initDevice() < 0 ) {
        exit(1);
      }
      
      pthread_create(&capture_thread, NULL, audioCapture, (void*)this);
      // Start recording thread
    }
    ~audioInput() {
      snd_pcm_close (pcm_handle);
    }

    void quitNow() {
      quit = true;
//      pthread_kill_other_threads_np();
      snd_pcm_close (pcm_handle);
    }
    
    union byte {
      unsigned char uchar_val;
      char char_val;
    };
    
    int index( int i ) {
      while( i < 0 ) {
        i += data_size;
      }
      return i % data_size;
    }
    
    float v(char b1, char b2) {
      byte b;
      b.char_val = b1;
      return (b2*256) + b.uchar_val;
    }
    static void* audioCapture(void* a) {
//      fprintf(stderr, "Capture thread started..\n");
      audioInput* ai = (audioInput*) a;
      
      float inv256 = 1.0 / 256.0;
      float inv256_2 = inv256*inv256;
      
      while( ! ai->quit ) {
        int n;
        if( ! ai->pause ) {
          while((n = snd_pcm_readi(ai->pcm_handle, ai->data, ai->read_frames)) < 0 ) {
            snd_pcm_prepare(ai->pcm_handle);
            fprintf(stderr, "Dropped audio data.\n");
          }
          
          byte b;
          int write_ptr, read_ptr;
          for( int i = 0; i < n; i++ ) {
            read_ptr = i * 2;
            write_ptr = ai->index(ai->data_write + i);
            b.char_val = ai->data[read_ptr];
            ai->float_data[write_ptr] = (float)ai->data[read_ptr+1]*inv256 + (float)b.uchar_val*inv256_2;
          }
          
          ai->data_end   = ai->data_write;
          ai->data_write = ai->index(ai->data_write+n);
          ai->data_start = ai->data_end - (ai->data_size - ai->write_padding);
          timeval t;
          gettimeofday(&t, NULL);
          ai->avg_write_interval = ai->avg_write_interval * 0.7 + timeDiff(ai->last_write, t) * 0.3;
          ai->last_write = t;
        }
        else {
          usleep(10000);
        }
      }
      fprintf(stderr, "Capture thread exiting.\n");
    }
    
    void getTimeSpan(float* time, int* index, int num) {
      float sh;
      if( pause ) {
        sh = 0;
      }
      else {
        timeval t;
        gettimeofday(&t, NULL);
        sh = timeDiff(last_write, t);
      }
      
//      float min_time = front_padding*dt > sh ? front_padding*dt : sh;
      float min_time = 0;
      float max_time = (data_end - data_start) * dt;
      
      for( int i=0; i<num; i++ ) {
        if( time[i] < min_time )
          time[i] = min_time;
        if( time[i] > max_time )
          time[i] = max_time;
        
/*        index[i] = data_end - (int)((time[i]-sh)*(float)rate);
        time[i] = sh + ((data_end - index[i]) * dt);*/
        index[i] = data_end - front_padding - (int)((time[i]-sh)*(float)rate);
        time[i] = sh + ((data_end - index[i] - front_padding) * dt);
      }
    }
    
    
    int initDevice() {
  
      stream = SND_PCM_STREAM_CAPTURE;
      /* Init pcm_name. Of course, later you */
      /* will make this configurable ;-)     */
      pcm_name = strdup("plughw:0,0");
      //pcm_name = strdup("hw:0,0");
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
  
      /* Set sample rate. If the exact rate is not supported */
      /* by the hardware, use nearest possible rate.         */ 
      exact_rate = rate;
      if (snd_pcm_hw_params_set_rate_near(pcm_handle, hwparams, (uint*)&exact_rate, 0) < 0) {  // note skipping the (uint*) here caused assertion failure in snd setup later!
        fprintf(stderr, "Error setting rate.\n");
        return(-1);
      }
      if (rate != exact_rate) {
        fprintf(stderr, "The rate %d Hz is not supported by your hardware.\n \
                        ==> Using %d Hz instead.\n", rate, exact_rate);
      }
  
      /* Set number of channels - note it's 1 otherwise later segfaults */
      if (snd_pcm_hw_params_set_channels(pcm_handle, hwparams, 1) < 0) {
        fprintf(stderr, "Error setting channels.\n");
        return(-1);
      }
  
      /* Set number of periods. Periods used to be called fragments. */ 
      if (snd_pcm_hw_params_set_periods(pcm_handle, hwparams, periods, 0) < 0) {
        fprintf(stderr, "Error setting periods.\n");
        return(-1);
      }
      /* Set buffer size (in frames). The resulting latency is given by */
      /* latency = periodsize * periods / (rate * bytes_per_frame)     */
      if (snd_pcm_hw_params_set_buffer_size_near(pcm_handle, hwparams, &exact_size) < 0) {
        fprintf(stderr, "Error setting buffersize.\n");
	return(-1);
      }
      if( size != exact_size ) {
        fprintf(stderr, "Buffer size %d is not supported, using %d instead.\n", (int)size, (int)exact_size);
      }

        /* Apply HW parameter settings to */
      /* PCM device and prepare device  */
      if (snd_pcm_hw_params(pcm_handle, hwparams) < 0) { // fails if earlier exact_rate not of uint* type, with uniformative error! ahb 5/16/15
        fprintf(stderr, "Error setting HW params.\n");
        return(-1);
	}
      return 1;
    }

};



class scene
{
  public:
    float dim[2];
    float aspect;
    Vector mouse;

    float side[4]; // t, b, r, l;
    float target_side[4];
    bool pause;
    float trigger[2];
    float trigger_width;
    int t_dir;
    audioInput* ai;
    float times[3];  // start, stop, trigger
    int indexes[3];
    GLuint tex;
    float line_width;        // Barnett added
  bool antialias;   // ditto

    scene() {
      t_dir = 1;
      side[0] = 0.3;
      side[1] = -0.3;
      side[2] = 0;
      side[3] = -0.1;
      for( int i=0; i<4; i++ )
        target_side[i] = side[i];
      
      trigger[0] = 0.0;  // y-value of trigger
      trigger[1] = -0.0002;  // x-value of trigger, was -0.05, Barnett
      trigger_width = 0.0002;     // what's this for? Barnett
      line_width = 2.0;          // start value
      antialias = true;          // Barnett default
      pause = false;
      ai = new audioInput();
      glGenTextures(1, &tex);
    }
    ~scene() {
    }
    
    float smooth(float a, float b, float s) {
      return a + (b-a)*s;
    }
  
    void setProjection() {
      float s = 0.3;
      glLoadIdentity();
      glOrtho(side[3], side[2], side[1], side[0], -10, 10);  // l r b t n f
    }
  
    void init(int argc, char** argv) {
      float white[] = {1.0, 1.0, 1.0, 1.0};
      float gray2[] = {0.25, 0.25, 0.25, 1.0};
    
      glLightfv(GL_LIGHT0, GL_SPECULAR, white);
      glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
      glLightfv(GL_LIGHT0, GL_AMBIENT, gray2);
      
      glLightfv(GL_LIGHT1, GL_SPECULAR, white);
      glLightfv(GL_LIGHT1, GL_DIFFUSE, white);
      glLightfv(GL_LIGHT1, GL_AMBIENT, gray2);
      
      GLfloat pos[4] = {0, 0, 20, 1};
      glLightfv (GL_LIGHT0, GL_POSITION, pos);
      GLfloat pos1[4] = {0, -10, -20, 1};
      glLightfv (GL_LIGHT1, GL_POSITION, pos1);
      
      glClearColor (0.0, 0.0, 0.0, 0.0);
    
      glMatrixMode (GL_MODELVIEW);
      glLoadIdentity ();
    }

    void printData() {
      
      times[0] = ai->data_size * ai->dt;
      times[1] = 0; // -side[2];
      ai->getTimeSpan(times, indexes, 3);
      for( int i=indexes[0]; i<indexes[1]; i++ ) {
        printf( "%f, %f\n", (i-indexes[0])*ai->dt, ai->float_data[ai->index(i)] );
      }
      printf("\n");
      fprintf(stderr, "Finished printing data.\n");
    }
  
    void drawPlot(float* data) {
      float dt = ai->dt;
      
//      if( ! pause ) {
        times[0] = -side[3];
        times[1] = -side[2];
        times[2] = -trigger[1];
        ai->getTimeSpan(times, indexes, 3);
        
        if( t_dir != 0 ) {
          int dx = 0;
          int tw = (int)(trigger_width / dt);
          if( t_dir < 0 ) {
            while( dx < 1000 && ! (ai->float_data[ai->index(indexes[2]+dx)] > trigger[0] && ai->float_data[ai->index(indexes[2]+dx+tw)] < trigger[0] )  ) {
              dx++;
            }
          }
          else {
            while( dx < 1000 && ! (ai->float_data[ai->index(indexes[2]+dx)] < trigger[0] && ai->float_data[ai->index(indexes[2]+dx+tw)] > trigger[0] )  ) {
              dx++;
            }
          }
          if( dx == 1000 )
            dx = 0;
          indexes[0] += dx;
          indexes[0] --;
          times[0] += dt;
          indexes[1] += dx;
          indexes[1]++;
        }
//      }
      
      
      
      
      int start = indexes[0];
      int stop = indexes[1];
      float startx = times[0];
      
      glLineWidth((int)line_width);    // Barnett
      if (antialias) {
	glEnable(GL_LINE_SMOOTH);    // Barnett anti-aliasing just for graph (line strip)
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      } else {
	glDisable(GL_BLEND);       // Barnett no anti-aliasing
	glDisable(GL_LINE_SMOOTH);
      }
      glDisable(GL_DEPTH_TEST);
      glDisable(GL_LIGHTING);

      glColor4f(0.4, 1.0, 0.4, 1);    // Barnett made brighter; was 0 0.6 0
      glBegin(GL_LINE_STRIP);
      float x = -startx;
      for( int i=start; i<=stop; i++ ) {
//        printf("i: %d   stop: %d\n", i, stop_index);
        glVertex2f(x, ai->float_data[ai->index(i)]);
        x += dt;
      }
      glEnd();
 
     
      glDisable(GL_LINE_SMOOTH);
      glBlendFunc (GL_SRC_ALPHA, GL_DST_ALPHA);  // these 5 lines moved from above line strip code
      glEnable (GL_BLEND);
      glEnable (GL_ALPHA_TEST);
      glAlphaFunc (GL_ALWAYS, 0.5);
      glEnable(GL_ALPHA);

     float pix_density = dim[0]/(side[2]-side[3]);
      if( ai->rate < pix_density ) {
        float c = (1-ai->rate/pix_density);
        glColor4f(c, c, 0, 1);
        glPointSize(2);
        glBegin(GL_POINTS);
        float x = -startx;
        for( int i=start; i<=stop; i++ ) {
          glVertex2f(x, ai->float_data[ai->index(i)]);
          x += dt;
        }
        glEnd();
      }
      glLineWidth(1);    // Barnett
      
      glColor4f(1, 1, 0, 0.3);
      glBegin(GL_LINES);
        glVertex2f(side[3], trigger[0]);
        glVertex2f(side[2], trigger[0]);
        glVertex2f(trigger[1], side[0]);
        glVertex2f(trigger[1], side[1]);
        glVertex2f(trigger[1]+trigger_width, side[0]);
        glVertex2f(trigger[1]+trigger_width, side[1]);
      glEnd();
      
      glColor4f(0.5, 0.5, 1, 0.5);
      glBegin(GL_LINES);
        glVertex2f(0, side[0]);
        glVertex2f(0, side[1]);
      glEnd();
    }

    void drawPlot() {
      drawPlot(ai->float_data);
      for( int i=0; i<4; i++ )
        side[i] = smooth(side[i], target_side[i], 0.3);
    }
    
    void scale(int ax, float s) {
      int s1 = ax*2;
      int s2 = s1+1;
      float w = target_side[s1]-target_side[s2];
      float dw = (w - w*s) / 2.0;
      target_side[s1] += dw;
      target_side[s2] -= dw;
    }
  void move(int ax, float x) {               // ax = 0 is y-direction,  ax = 1 is x-direction (weird)
      int s1 = ax*2;
      int s2 = s1+1;
      float w = target_side[s2]-target_side[s1];
      float dx = x*w;
      target_side[s1] += dx;
      target_side[s2] += dx;
    }
};
static scene scn;


void drawText(float x, float y, char *string) {
  int len, i;
  glRasterPos2f(x, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++)
  {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
  }
}


void drawGrid() {

  char* str = new char[50];
  char* fmt = new char[10];
  glLineWidth(1.0);
  glBlendFunc (GL_SRC_ALPHA, GL_DST_ALPHA);
  glEnable (GL_BLEND);
  glEnable (GL_ALPHA_TEST);
  glAlphaFunc (GL_ALWAYS, 0.5);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_ALPHA);
  glDisable(GL_LIGHTING);
  
  int d = -1;
  float dt = scn.side[2] - scn.side[3];
  for( float ivl = 10; ivl > .000001; ivl *= 0.1 ) {
    if( ivl < dt*5 && ivl > dt/100 ) {
      int min = (int)(scn.side[3] / ivl)-1;
      int max = (int)(scn.side[2] / ivl)+1;
      max = max > 0 ? 0 : max;
      float alpha = 5*ivl/dt;
      glColor4f(1, 1, 1, alpha > 0.1 ? 0.1 : alpha);
      glBegin(GL_LINES);
      for( int i=min; i<max; i++ ) {
        glVertex2f(i*ivl, scn.side[0] < 1 ? scn.side[0] : 1);
        glVertex2f(i*ivl, scn.side[1] > -1 ? scn.side[1] : -1);
      }
      glEnd();
      if( ivl > dt/20.0 ) {
        glColor4f(1, 1, 1, alpha > 0.5 ? 0.4 : alpha - 0.1);
        for( float i=min; i<max; i++ ) {
          sprintf(str, "%g", -i*ivl);
          drawText(i*ivl, scn.side[1] * 0.99 + scn.side[0] * 0.01, str);
        }
      }
    }
    d++;
  }

  float dv = scn.side[0] - scn.side[1];
  for( float ivl = 10.0; ivl > .000001; ivl *= 0.1 ) {
    if( ivl < dv*5.0 && ivl > dv/200.0 ) {
      
      int min = (int)(scn.side[1] / ivl)-1;
      int max = (int)(scn.side[0] / ivl)+1;
      min = scn.side[1] < -1 ? (-1/ivl) : min;
      max = scn.side[0] > 1 ? (1/ivl) : max;
      
      float alpha = 5*ivl/dv;
      glColor4f(1, 1, 1, alpha > 0.1 ? 0.1 : alpha);
      glBegin(GL_LINES);
      for( int i=min; i<=max; i++ ) {
        glVertex2f(scn.side[2] > 0 ? 0 : scn.side[2], i*ivl);
        glVertex2f(scn.side[3], i*ivl);
      }
      glEnd();
      if( ivl > dv/20.0 ) {
        glColor4f(1, 1, 1, alpha > 0.5 ? 0.4 : alpha - 0.1);
        for( float i=min; i<=max; i++ ) {
          sprintf(str, "%g", i*ivl);
          drawText(scn.side[3] * 0.99 + scn.side[2] * 0.01, i*ivl, str);
        }
      }
    }
  }
  
  
}


void drawUI() {
  glDisable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
      glLoadIdentity();
      glOrtho(-1, 1, -1, 1, -1000, 1000);
      glColor4f(1, 1, 1, 0.2);
      glBegin(GL_LINES);
        glVertex2f(0, 0.1);
        glVertex2f(0, -0.1);
        glVertex2f(0.1, 0);
        glVertex2f(-0.1, 0);
      glEnd();
      glColor4f(1, 1, 1, 1);
//      drawText(-0.95, 0.95, "UI");
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}


  
void display()
{

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  scn.setProjection();
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
    glLoadIdentity();
    drawGrid();
    scn.drawPlot();
  glPopMatrix();
  
  drawUI();    
  glFinish();
  glutSwapBuffers();
}

void idle(void)
{
  glutPostRedisplay();
}

void special(int key, int xPos, int yPos) {
  if( key == 101 ) { // up
    scn.move(0, -0.2);
  }
  else if( key == 103 ) { // dn       - Barnett reversed signs for more natural interface
    scn.move(0, 0.2);
  }
  else if( key == 100 ) { // lt
    scn.move(1, 0.2);
  }
  else if( key == 102 ) { // rt
    scn.move(1, -0.2);
  }
  else if( key == 104 ) { // pgup
    scn.scale(0, 1.1);
  }
  else if( key == 105 ) { // pgdn
    scn.scale(0, 1/1.1);
  }
  else if( key == 106 ) { // home
    scn.scale(1, 1.1);
  }
  else if( key == 107 ) { // end
    scn.scale(1, 1/1.1);
  }
  else {
    fprintf(stderr, "pressed special key %d\n", key);
  }
}

void keyboard(unsigned char key, int xPos, int yPos)
{
  if( key == 27 || key == 'q') {  // esc or q to quit
    scn.ai->quitNow();
    exit(0);
  }
  else if( key == ' ' ) {
    scn.pause = ! scn.pause;
    scn.ai->pause = ! scn.ai->pause;
  }
  else if( key == 't' ) {
    scn.t_dir = ((scn.t_dir+2) % 3) - 1;
  }
  else if( key == 'p' ) {
    scn.printData();
  }
  else if( key == '2' ) { // Barnett line width
    scn.line_width = min(10.0, scn.line_width+1.0);
  }
  else if( key == '1' ) { // Barnett
    scn.line_width = max(1.0, scn.line_width-1.0);
  }
  else if( key == 'l' ) { // Barnett antialiasing toggle
    scn.antialias = !scn.antialias;
  }
  else if( key == 9 ) { // Tab         - Barnett: to move in x to align right side with trigger
    scn.move(1, -(scn.side[2] - scn.trigger[1])/(scn.side[3]-scn.side[2])); 
  }
  else {
    fprintf(stderr, "pressed key %d\n", (int)key);
  }
}

void reshape(int w, int h)
{
  scn.dim[0] = w;
  scn.dim[1] = h;
//  fprintf(stderr, "Setting w=%d, h=%d\n", w, h);
  glViewport(0, 0, w, h);
}

void mouse(int button, int state, int x, int y)
{
  scn.mouse[0] = x;
  scn.mouse[1] = y;
  scn.mouse[2] = button;
  if( button == 3 && state == 1 ) {           // I don't get this - Barnett
    scn.scale(0, 1.1);
    scn.scale(1, 1.1);
  }
  else if( button == 4 && state == 1 ) {
    scn.scale(0, 1/1.1);
    scn.scale(1, 1/1.1);
  }
  else {
//    printf("Mouse event: button: %d  state: %d  pos: %d, %d\n", button, state, x, y);
  }
}

void motion(int x, int y)
{
  int dx = (int)(x-scn.mouse[0]);
  int dy = (int)(y-scn.mouse[1]);
  if( scn.mouse[2] == 0 ) {
    scn.move(0, -((float)dy)/scn.dim[1]);
    scn.move(1, ((float)dx)/scn.dim[0]);
  }
  else if( scn.mouse[2] == 2 ) {
    scn.scale(0, 1.0 - dy/200.0);  // made 4x less sensitive, Barnett
    scn.scale(1, 1.0 + dx/200.0);
  }
  else if( scn.mouse[2] == 1 ) {
    scn.trigger[0] -= dy/scn.dim[1] * (scn.side[0]-scn.side[1]);
    scn.trigger[1] += dx/scn.dim[0] * (scn.side[2]-scn.side[3]);
  }
  scn.mouse[0] = x;
  scn.mouse[1] = y;
}



int main(int argc, char** argv)
{
  printf("glscope by Luke Campagnola 2005, tweaked by Alex Barnett 5/16/15\n");
  printf("\nPgUp/PgDn vertical zoom; mousewheel overall zoom; arrows pan\n");
  printf("Home/End horiz zoom; Space pause; 1/2 line width; l antialias\n");
  printf("p data to stdout; q quit\n");
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(800, 600);
//  glutInitWindowPosition(0, 0);
  int mainWindow = glutCreateWindow("glScope: by Luke Campagnola, modified by Alex Barnett");
//  glutFullScreen();
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutSpecialFunc(special);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutIdleFunc(idle);
  scn.init(argc, argv);
  glutMainLoop();
  return 0;
}

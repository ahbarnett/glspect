# glspect

Real-time audio spectrogram in OpenGL

Alex Barnett	 	     		Mathematics Dept, Dartmouth College
	 	   	  		December 2010
					(updated Jan 2015. github 3/15/17)
					ahb@math.dartmouth.edu

This is a simple real-time spectrogram using OpenGL, GLUT, and ALSA.
It displays a scrolling window showing a "player piano roll" of the recent audio received by the microphone, with frequency on the vertical axis, and time on the horizontal.
It is excellent for demonstrating many acoustic phenomena such as vocal formants, resonance, doppler shifts, vibrational frequencies of instruments and metallic objects, timbre (harmonic content), etc.
The color range and contrast, as well as the time scroll rate, are adjustable
in real time. The frequency range is command-line adjustable in powers of 2.
Smaller windows at the bottom display real-time graphs of the spectrum and the
signal.

### Installation

This assumes a linux system.
Ensure you have the dependencies installed:


Run the Makefile via "make". You will need the math
library, single-precision FFTW libraries (RPM packages fftw and
fftw-devel), and ALSA (alsa-lib, alsa-lib-devel).  You will also want
to enable vSync (vertical refresh) in your graphics card settings (for
my machine, this menu is reached via nvidia-settings).  If your
refresh rate is not 60 Hz, you should change this manually in the code
(sorry; this will be made into an option in future).

(Jan 2015: I just discovered glDrawPixels is deprecated in OpenGL >3.0,
which is bad news - if someone figures out how to replace it, let me know).


### Usage:

```bash
glspect  [-f] [-v] [-sf <scroll_factor>] [-w <windowtype>] [-t twowinsize]

Command line arguments:
windowtype =    0 (no window) (will be crappy)
                1 (Hann)
                2 (Gaussian trunc at +-4sigma) (default, recommended)
scroll_factor = 1,2,... #. How many vSyncs (@ 60Hz) to wait per scroll pixel
                  (default 1)
twowinsize = 11,12,...,16 is the power of 2 giving FFT win_size N (default 13)
        (Note: this controls the vertical frequency resolution and range)

Keys & mouse:   arrows or middle button drag - brightness/contrast
       		i - step through colormaps (B/W, inverse B/W, color)
                q or Esc - quit
                [ and ] - control horizontal scroll factor (rate)
```

### Bug reports

Please submit an Issue to github, or contact me at the above email address.


### To do

* marking and playback of selections, saving of recent audio  
* real-time change of frequency range (win_size) ... but annoying  
* better docs
* ports to Mac OSX, Windows  

Please submit a pull request if you fix one of these!


### License, provenance, notes

With default scroll_factor=2, CPU usage is less than 50% of one core of a modern intel i7 CPU.
With scroll_factor=1 (fastest scrolling), this may go up to 75%. If your CPU or GPU is not as fast as this, you might want to shrink n_f and/or n_tw in the code.

The code is somewhat based upon that of `glScope` by Luke Campagnola (2005).
It is distributed under a completely free license; this means you can do
absolutely anything you want with this code. I'd appreciate if you
credit me where appropriate, though.

It is also influenced by `baudline` (an amazing tool, which unfortunately has the time and frequency axes flipped).

Here's some other info taken from Luke's `glScope README` file:

```
Since there is no interface to change the audio input device that is
displayed, you'll actually have to change some code if you want any
other channel besides the default recording device. The important line
in the code looks like this: 

pcm_name = strdup("plughw:0,0");

I guess you'll have to read about ALSA to figure out how to change 
that string.
```

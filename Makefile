# linux makefile for glspect.
# Barnett 12/20/10. Tidied 1/25/15. Tweaked for github 3/15/17.

LD_FLAGS =  
LD_LIBS = -lasound -lfftw3f -lm -pthread  # note single-prec FFTW library
CXX = g++
CXX_FLAGS = -O2

all:	glspect

glspect: glspect.cpp Makefile
	${CXX} ${CXX_FLAGS} glspect.cpp -lglut -lGL ${LD_LIBS} -o glspect

clean:
	rm -f core *.o glspect

tar:
	(cd ..; tar cvfz glspect.tgz glspect/README.md glspect/Makefile glspect/glspect.cpp glspect/glspect glspect/screenshot.png)

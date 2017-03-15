# glSpect, etc, Makefile, taken from glscope Makefile
# Barnett 12/20/10. Tidied 1/25/15

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
	(cd ..; tar cvfz glSpect.tar.gz glSpect/README glSpect/Makefile glSpect/glspect.cpp glSpect/glspect)

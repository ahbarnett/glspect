# luke canpagnola 2005, tweaked alex barnett 5/16/15
LD_FLAGS =  
LD_LIBS = -lasound -lm -pthread
CXX = g++
#CXX_FLAGS = -g -Wno-deprecated
CXX_FLAGS = -O2

all: glscope

glscope: glscope.cpp Makefile
	${CXX} ${CXX_FLAGS} glscope.cpp -lglut -lGL ${LD_LIBS} -o glscope

clean:
	rm -f core *.o glscope

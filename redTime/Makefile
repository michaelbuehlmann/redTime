# CCFLAGS=-O3 -fopenmp -std=c++11 -g
CCFLAGS=-O0 -std=c++11 -g
PATHS=-I/usr/include -L/usr/lib/x86_64-linux-gnu/
LIBS=-lgsl -lgslcblas

# modern camb?
CCFLAGS += -DMODERN_CAMB

redTime: redTime.cc AU_cosmological_parameters.h AU_interp.h AU_tabfun.h
	g++ redTime.cc -o redTime $(CCFLAGS) $(PATHS) $(LIBS)

clean:
	$(RM) redTime



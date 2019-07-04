CCFLAGS=-O3 -fopenmp
PATHS=-I/usr/include -L/usr/lib/x86_64-linux-gnu/
LIBS=-lgsl -lgslcblas

redTime: redTime.cc AU_cosmological_parameters.h AU_interp.h AU_tabfun.h
	g++ redTime.cc -o redTime $(CCFLAGS) $(PATHS) $(LIBS)

clean:
	$(RM) redTime



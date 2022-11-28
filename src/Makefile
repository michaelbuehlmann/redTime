CFLAGS=-O3 -g -Wall -Wpedantic
CCFLAGS=-O3 -fopenmp -std=c++11 -g
# CCFLAGS=-O0 -std=c++11 -g
PATHS=-I/usr/include -L/usr/lib/x86_64-linux-gnu/
LIBS=-lgsl -lgslcblas -lm

# modern camb?
# CCFLAGS += -DMODERN_CAMB

.PHONY: all
all: redTime convertPk

redTime: redTime.cc AU_cosmological_parameters.h AU_interp.h AU_tabfun.h
	g++ redTime.cc -o redTime $(CCFLAGS) $(PATHS) $(LIBS)

convertPkFull: convert_pk.c
	gcc convert_pk.c -o convertPkFull $(CFLAGS) $(PATHS) $(LIBS)

convertPt: convert_pt.c
	gcc convert_pt.c -o convertPT $(CFLAGS) $(PATHS) $(LIBS)

clean:
	$(RM) redTime
	$(RM) convertPk



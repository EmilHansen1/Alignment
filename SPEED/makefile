CC = gcc
CFLAGS = $(shell gsl-config --cflags) -fopenmp -Wall -Ofast -std=gnu11
LDLIBS = $(shell gsl-config --libs) -lm

default: main
	./$<

main: main.o units.h planar_rotor.o utils.o
utils.o: utils.h units.h
planar_rotor.o: planar_rotor.h units.h utils.o

clean:
	$(RM) main *.o *.out
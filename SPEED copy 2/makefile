CC = gcc
CFLAGS = $(shell gsl-config --cflags) -fopenmp -Wall -Ofast -std=gnu11
LDLIBS = $(shell gsl-config --libs) -lm -lgomp

default: main
	./$<
	python3 plot.py

main: main.o units.o planar_rotor.o utils.o input_loader.o #Linear_rotor.o
units.o: units.h
utils.o: utils.h units.o
planar_rotor.o: planar_rotor.h units.o utils.o
input_loader.o: input_loader.h planar_rotor.o units.o utils.o
#linear_rotor.o: linear_rotor.h units.h utils.o

clean:
	$(RM) main *.o *.out
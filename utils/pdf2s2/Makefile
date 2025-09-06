LIB = -lm -lgsl -lgslcblas
CC = gcc
DEBUG = -g

all : lattice.o main.o
	$(CC) -o pdf2s2_v2 *.o $(LIB) $(DEBUG)

lattice.o : lattice.c lattice.h
	gcc -c lattice.c $(DEBUG)

main.o : main.c
	gcc -c main.c $(DEBUG)

clean :
	rm *.o

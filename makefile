# Makefile for wobble

eigenflag = -I ~/include/eigen-3.4.0

all: wobble

wobble: src/wobble.o
	gcc src/wobble.o -o wobble -lstdc++ -lm

src/wobble.o: src/wobble.cpp
	gcc -g $(eigenflag) -c src/wobble.cpp -o src/wobble.o 

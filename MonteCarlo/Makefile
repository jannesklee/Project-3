CC = c++
DEBUG = -g
CFLAGS = -Wall -larmadillo -llapack -lblas -fopenmp $(DEBUG)
DEPS = singleparticle.h manybody.h lib.h
OBJ = main.o singleparticle.o manybody.o lib.o

%.o: $(DEPS)
	$(CC) $(CFLAGS)-c -o $@ $<

main.out: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ 

all: main.out 

.PHONY: clean

clean:
	$(RM) *.o *.out


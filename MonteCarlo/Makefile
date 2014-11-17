CC = c++
CFLAGS = -I -Wall -larmadillo -llapack -lblas -fopenmp
DEPS = singleparticle.h manybody.h lib.h
OBJ = main.o singleparticle.o manybody.o lib.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main.out: $(OBJ)
	c++ -o $@ $^ $(CFLAGS)

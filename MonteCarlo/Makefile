CC = c++
DEBUG = -g
CFLAGS = -Wall -larmadillo -llapack -lblas -fopenmp -O2 $(DEBUG)
DEPS = singleparticle.h manybody.h lib.h
OBJ = main.o singleparticle.o manybody.o lib.o

%.o: %.cpp $(DEPS)
	LC_MESSAGES=C $(CC) -c $(CFLAGS) -o $@ $<

main.out: $(OBJ)
	LC_MESSAGES=C $(CC) $(CFLAGS) -o $@ $^ 

all: main.out 

.PHONY: clean

clean:
	$(RM) *.o *.out


CC = c++
CFLAGS = -I -Wall -larmadillo -llapack -lblas
DEPS = singleparticle.h manybody.h
OBJ = main.o singleparticle.o manybody.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main.out: $(OBJ)
	c++ -o $@ $^ $(CFLAGS)

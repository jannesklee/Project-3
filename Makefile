TARGET = main.out
CC = c++
CFLAGS = -Wall -Wextra -Isrc -larmadillo -llapack -lblas -fopenmp -O2 $(DEBUG)
DEBUG = -g

SRCDIR = src
OBJDIR = obj
#BINDIR = bin

SOURCES = $(wildcard $(SRCDIR)/**/*.cpp $(SRCDIR)/*.cpp)
DEPS = $(wildcard $(SRCDIR)/**/*.h $(SRCDIR)/*.h)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o,$(SOURCES))

# not used yet
#TEST_SRC = $(wildcard tests/*_tests.cpp)
#TESTS = $(patsubst %.cpp,%,$(TEST_SRC))

# compiling
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	LC_MESSAGES=C $(CC) -c $(CFLAGS) -o $@ $<

# linking
$(TARGET): $(OBJECTS)
	LC_MESSAGES=C $(CC) $(CFLAGS) -o $@ $^ 

all: main.out 

# cleaning procedure
.PHONY: clean
clean:
	$(RM) $(OBJECTS) $(TARGET)


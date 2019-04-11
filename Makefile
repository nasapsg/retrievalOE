# the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
CFLAGS = -Wall -Os -pedantic
LIBS = -lm

# Output directory
PROGS = oem

# Name of the executables:
all: $(PROGS)

oem: oemtest.c oem.c oem.h 
	$(CC) $(CFLAGS) -lm -o $@ oemtest.c oem.c

clean:
	rm -f $(PROGS)

CC = g++
CFLAGS = -fopenmp -O3 -g

all: tcrmatch

tcrmatch: src/tcrmatch.cpp
	$(CC) $(CFLAGS) -o tcrmatch src/tcrmatch.cpp

clean:
	rm tcrmatch

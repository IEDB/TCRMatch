CC = g++
CFLAGS = -std=c++11 -fopenmp -O3 -g

all: tcrmatch

tcrmatch: src/tcrmatch.cpp
	$(CC) $(CFLAGS) -o tcrmatch src/tcrmatch.cpp

clean:
	rm tcrmatch

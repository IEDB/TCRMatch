CC = g++
CFLAGS = -std=c++11 -fopenmp -O3 -g

all: tcrmatch

tcrmatch: src/main.cpp src/tcrmatch.cpp
	$(CC) $(CFLAGS) -o tcrmatch src/main.cpp 

test: tcrmatch
	bash tests/run_tests.sh

clean:
	rm tcrmatch

.PHONY: all test clean


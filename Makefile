CC = g++
CFLAGS = -std=c++11 -fopenmp -O3 -g

all: tcrmatch

tcrmatch: src/*.cpp
	$(CC) $(CFLAGS) -o tcrmatch src/*.cpp 

clean:
	rm tcrmatch


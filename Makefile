CC = g++
CFLAGS = -fopenmp -O3 -g

all: tcrmatch

tcrmatch: src/main.cpp
	$(CC) $(CFLAGS) -o tcrmatch src/main.cpp

clean:
	rm tcrmatch

CC = g++
CFLAGS = -fopenmp -O3 -g

all: tcrmatch

tcrmatch: src/main.cpp
	$(CC) $(CFLAGS) -o tcrmatch src/main.cpp
	$(CC) $(CFLAGS) -o test tests/test.cpp tests/catch_amalgamated.cpp

clean:
	rm tcrmatch
	rm test

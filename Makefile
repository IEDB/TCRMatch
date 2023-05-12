CC = g++
CFLAGS = -std=c++11 -fopenmp -O3 -g

all: tcrmatch

tcrmatch: src/main.cpp src/tcrmatch.cpp
	$(CC) $(CFLAGS) -o tcrmatch src/main.cpp

test: tests/test.cpp tests/catch_amalgamated.cpp src/main.cpp src/tcrmatch.cpp
	$(CC) $(CFLAGS) -o test tests/test.cpp tests/catch_amalgamated.cpp src/main.cpp src/tcrmatch.cpp

clean:
	rm tcrmatch
	rm test

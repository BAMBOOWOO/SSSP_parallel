CC=mpic++
CFLAGS=-O2 -fopenmp -std=c++14 -I./lib/usr/include

main: main.o common.o parallel.o
	$(CC) $(CFLAGS) main.o common.o parallel.o -o main

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

common.o: common.cpp
	$(CC) $(CFLAGS) -c common.cpp

parallel.o: parallel.cpp
	$(CC) $(CFLAGS) -c parallel.cpp
clean:
	rm -rf *.o main
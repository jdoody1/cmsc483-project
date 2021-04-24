CC=gcc
CFLAGS=-g -fopenmp
main: grauer_translation.o
	$(CC) -o grauer_translation.o -lm

parallel:
	$(CC) $(CFLAGS) grauer_parallel.c -o grauer_parallel.o -lm

clean:
	rm *.o main

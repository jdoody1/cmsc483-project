CC=gcc
CFLAGS=-g -fopenmp
main: grauer_translation.o
	$(CC) grauer_translation.c -o grauer_translation.o -lm

parallel:
	$(CC) $(CFLAGS) grauer_parallel.c -o grauer_parallel.o -lm

clean:
	rm *.o main

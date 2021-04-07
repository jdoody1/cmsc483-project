CC=gcc

main: grauer_translation.o
	$(CC) -o main grauer_translation.o -lm

clean:
	rm *.o main
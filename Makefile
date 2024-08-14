CC=clang
CFLAGS=-Wall -pedantic

main: main.c lib.o
	$(CC) $(CFLAGS) -o main main.c lib.o

lib.o: lib.c
	$(CC) $(CFLAGS) -o lib.o -c lib.c
	

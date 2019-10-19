all: nbody_sh1 sphere2 find_binaries2 find_escapees

# Which compiler
#CC = gcc
CC = g++

# Options for development
#CFLAGS = -g -Wall

# Options for release
CFLAGS = -O -Wall

nbody_sh1: nbody_sh1.o
	$(CC) -o nbody_sh1 nbody_sh1.o

nbody_sh1.o: nbody_sh1.c
	$(CC) $(CFLAGS) -c nbody_sh1.c


sphere2: sphere2.o
	$(CC) -o sphere2 sphere2.o

sphere2.o: sphere2.c
	$(CC) $(CFLAGS) -c sphere2.c


find_binaries2: find_binaries2.o
	$(CC) -o find_binaries2 find_binaries2.o

find_binaries2.o: find_binaries2.c
	$(CC) $(CFLAGS) -c find_binaries2.c


find_escapees: find_escapees.o
	$(CC) -o find_escapees find_escapees.o

find_escapees.o: ../find_escapees/main.c
	$(CC) $(CFLAGS) -c ../find_escapees/main.c

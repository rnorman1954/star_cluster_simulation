all: nbody1 createCluster findBinaries findEscapees

# Which compiler
#CC = gcc
CC = g++

# Options for development
#CFLAGS = -g -Wall

# Options for release
CFLAGS = -O -Wall

#----------------------------------
# Compile the various files

nbody1: nbody1.o
	$(CC) -o build/nbody1 nbody1.o

nbody1.o: nbody1.cpp
	$(CC) $(CFLAGS) -c nbody1.cpp


createCluster: createCluster.o
	$(CC) -o build/createCluster createCluster.o

createCluster.o: createCluster.cpp
	$(CC) $(CFLAGS) -c createCluster.cpp


findBinaries: findBinaries.o
	$(CC) -o build/findBinaries findBinaries.o

findBinaries.o: findBinaries.cpp
	$(CC) $(CFLAGS) -c findBinaries.cpp


findEscapees: findEscapees.o
	$(CC) -o build/findEscapees findEscapees.o

findEscapees.o: findEscapees.cpp
	$(CC) $(CFLAGS) -c findEscapees.cpp


analyseMassDist: analyseMassDist.o
	$(CC) -o build/analyseMassDist analyseMassDist.o

analyseMassDist.o: analyseMassDist.cpp
	$(CC) $(CFLAGS) -c analyseMassDist.cpp

#=========================================
# Clean up all of the object files
.PHONY: clean
clean:
	rm -f *.o

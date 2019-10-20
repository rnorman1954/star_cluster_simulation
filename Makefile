all: nbody1 createCluster findBinaries findEscapees analyseMassDist createAnimationFiles convert

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

createAnimationFiles: createAnimationFiles.o
	$(CC) -o build/createAnimationFiles createAnimationFiles.o

createAnimationFiles.o: createAnimationFiles.cpp
	$(CC) $(CFLAGS) -c createAnimationFiles.cpp

convert: convert.o
	$(CC) -o build/convert convert.o

convert.o: convert.cpp
	$(CC) $(CFLAGS) -c convert.cpp

#=========================================
# Clean up all of the object files
.PHONY: clean
clean:
	rm -f *.o

#-*-makefile-*-

SEQAN_LIB = .
CXX=g++ -std=c++14
CC=$(CXX)

WARN= -W -Wall
FLAGS=-I$(SEQAN_LIB) -DSEQAN_HAS_ZLIB=1 -fopenmp

# Flags for optimization OR for debugging
RELEASE_FLAGS=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
DEBUG_FLAGS=-g -fno-inline -DDEBUG -DBOUNDS_CHECK -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# Replace RELEASE_FLAGS by DEBUG_FLAGS for debugging
CXXFLAGS = $(WARN) $(RELEASE_FLAGS) $(INCLUDE)
#CC = g++ -fno-merge-constants

LDLIBS=-lz -lrt -pthread -lgomp

.cc.o:LRcaller.cc LRcaller.cc
	$(CC) -c $(CXXFLAGS) $(FLAGS) $(TOOLS) $< -o $@


all: LRcaller

LRcaller: LRcaller.o $(OBJECTS)
	$(CXX) $(INCLUDES) $(OBJECTS) LRcaller.o $(LDFLAGS) $(LDLIBS) -o LRcaller

LRcaller.o: LRcaller.cc

clean:
	rm -f *.o LRcaller

default:
	all

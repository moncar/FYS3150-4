CC=g++
CFLAGS=-c -O3 -Wall -fopenmp -std=c++11
LDFLAGS=-lgomp
SOURCES=main.cpp lib.cpp outputFile.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl5.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

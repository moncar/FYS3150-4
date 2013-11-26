CC=g++
CFLAGS=-c -o3 -Wall -fopenmp -std=c++0x
LDFLAGS=-lgomp
SOURCES=main.cpp lib.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl5.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

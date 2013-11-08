CC=g++
CFLAGS=-c -Wall -O3 -std=gnu++11
LDFLAGS=-larmadillo
SOURCES=main.cpp outputFile.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl4.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

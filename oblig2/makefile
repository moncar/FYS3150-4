CC=g++
CFLAGS=-c -Wall -O3
LDFLAGS=-larmadillo
SOURCES=main.cpp jacobirot.cpp lib.cpp outputFile.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl2.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

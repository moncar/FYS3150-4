CC=g++
CFLAGS=-c -Wall -O1 -std=gnu++11
LDFLAGS=-lgsl -larmadillo
SOURCES=main.cpp outputFile.cpp PDESolvers.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl4.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

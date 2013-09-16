CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=o1-mainprog.cpp functions.cpp lumodule.cpp outputFile.cpp lib.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

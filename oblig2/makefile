CC=g++
CFLAGS=-c -Wall
LDFLAGS=-larmadillo
SOURCES=main.cpp jacobirot.cpp lib.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl2.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

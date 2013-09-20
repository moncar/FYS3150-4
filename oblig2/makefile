CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp jacobirot.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl2.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

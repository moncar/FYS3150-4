CC=g++
CFLAGS=-c -Wall -fopenmp
LDFLAGS=-lgomp
SOURCES=main.cpp lib.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl5.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

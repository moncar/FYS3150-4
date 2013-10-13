CC=g++
CFLAGS=-c -Wall -O3 -std=gnu++11
LDFLAGS=-larmadillo
SOURCES=main.cpp outputFile.cpp solver_EC.cpp solver_RK4.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl3.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

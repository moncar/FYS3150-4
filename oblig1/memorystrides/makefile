CC=g++
CFLAGS=-c -Wall -O3 -march=native 
LDFLAGS=-larmadillo
SOURCES=memorystrides.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=memstrides.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

CC=g++
CFLAGS=-c -O3 -Wall -fopenmp -std=c++11
LDFLAGS=-lgomp -lgsl -lgslcblas
SOURCES=main.cpp lib.cpp outputFile.cpp wavefuncs.cpp mainGaussLegendre.cpp mainGaussHermite.cpp mainMonteCarloBruteForce.cpp mainMonteCarloImportanceSampling.cpp mainMonteCarloVMC1.cpp extfuncGaussHermite.cpp miscfuncs.cpp
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=obl5.x

all: $(SOURCES) $(EXECUTABLE)
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

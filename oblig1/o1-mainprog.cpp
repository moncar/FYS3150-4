// 2013-09-16 FYS4150
// Marius Berge Eide

// Project 1
// Numerical differentiation

#include <iostream>
#include <cmath>
#include <stdlib.h>

#include "diffsolver.h"
//#include "errorchecker.h"

using namespace std;

// ASSIGN FUNCTIONS AND DATATYPES
// int 

// TO DO -> replace functions with INLINE FUNCTIONS
double fSourceFunc(double x) {
    return ( 100.0 * exp(-10. * x) );
}

double fSolutionFunc(double x) {
    return ( 1. - (1 - exp(-10.)) * x - exp(-10. * x) );
}

inline double fRelError(double fCorr, double fDeviation) {
    return (log10( abs( (fDeviation - fCorr)/fCorr ) ));
}



int main(int argc, char* argv[]) {
    // Define number of iterations
    int N = 0;
    if (!argv[1]) { N = 100; }
    else { N =  atoi(argv[1]); }
    // Set limits and step size
    const double fXmin = 0.0;
    const double fXmax = 1.0;
    const double fHstep= (fXmax - fXmin)/double(N);

    double *fX= new double[N];
    double *fEstX= new double[N];
    double *fCorrX= new double[N];
    double *fErr= new double[N];
    
    for (int iii=0; iii<N; iii++) {
        fX[iii]    = fXmin + iii*fHstep;
        fEstX[iii] = fSourceFunc(fX[iii]);      // to be sent to SOLVER
        fCorrX[iii]= fSolutionFunc(fX[iii]);
    }

    // Send array adress holding source function to SOLVER
    // which modifies the array
    solver(fEstX, fHstep, N);

    // Estimate relative error;
    // using results in array and comparing these values
    // with results from analytical expression evaluated at the same points
    for (int iii=0; iii<N; iii++) {
        fErr[iii] = fRelError(fCorrX[iii], fEstX[iii]);
        
        cout << fX[iii] << "\t" << fCorrX[iii] << "\t" << fEstX[iii] << "\t" << fErr[iii] << endl;
    }

    delete [] fX;       // housekeeping
    delete [] fEstX;
    delete [] fCorrX;
    delete [] fErr;
}


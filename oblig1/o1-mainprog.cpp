// Project 1
// Numerical differentiation
// Deadline Mon 16 September 12pm

#include <iostream>
#include <cmath>

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


int main() {
    const int    N     = 150;
    const double fXmin = 0.0;
    const double fXmax = 1.0;
    const double fHstep= (fXmax - fXmin)/double(N);
    cout << fHstep << endl;

    double *fX;
    double *fEstX;
    double *fCorrX;
    double *fErr;
    fX = new double[N];
    fEstX  = new double[N];
    fCorrX = new double[N];
    fErr   = new double[N];

    // Send empty array to solver, receive a filled one

    fX[0] = 0.12345; 
    solver(&fX, fHstep, fSourceFunc);
    cout << fX[0] << endl;
    // Estimate relative error;
    // using results in array and comparing these values
    // with results from analytical expression evaluated at the same points
    
    for (int iii=0; iii<N; iii++) {
        fErr[iii] = log10(abs( (fEstX[iii] - fCorrX[iii])/fCorrX[iii] ));
    }

    delete [] fX;       // housekeeping
    delete [] fEstX;
    delete [] fCorrX;
    delete [] fErr;
}


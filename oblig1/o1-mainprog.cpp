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
    const double fXmin = 0.001;
    const double fXmax = 0.999;
    const double fHstep= (fXmax - fXmin)/double(N);

    double *fX= new double[N];     // array: x-values
    double *fEstX= new double[N];   // array: source-func(x) -> SOLVER
    double *fCorrX= new double[N];  // array: answer(x)
    double *fErr= new double[N];    // array: to hold difference Est vs Corr
    
    for (int iii=0; iii<N; iii++) {
        fX[iii]    = fXmin + iii*fHstep; 
        fEstX[iii] = pow(fHstep,2) * fSourceFunc(fX[iii]); // h^2 * f(x)
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



    // PART 2:
    // LU DECOMPOSITIONING
    
    double *fEstX2= new double[N];  // array: source-func(x) -> LU decomp
   
    // Initialise array 
    for (int iii=0; iii<N; iii++) {
        fEstX2[iii] = pow(fHstep,2) * fSourceFunc(fX[iii]); // h^2 * f(x) 
    }

    // Apply LU decompositioning of differentiation tridiag matrix
    // which is calculated in LUmethod,
    // solve for x in "Ax = b" where A=LU and b=fEstX2
    // Returns x in place of fEstX2. 

    LUmethod(N, fEstX2);

    // Save result to file
    if (argc <= 2) { 
        cout << "ERROR! Missing output file for LU decomp result" << endl;
        exit(1);
    }
    else {
        outputFile(N, fX, fEstX2, &argv[2]);
    }
    
    // PART 3:
    // Matrix operations

    delete [] fX;       // housekeeping
    delete [] fEstX;
    delete [] fEstX2;
    delete [] fCorrX;
    delete [] fErr;



}


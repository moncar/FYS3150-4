// EULER CROMER SOLVER ROUTINE

// Swallows functions, initial values and solution arrays

//#include <iostream>
//#include <cmath>
#include <armadillo>

#include "headers.h"

using namespace std;
using namespace arma;

void solver_EC(   arma::vec (*vFuncIn)(double, arma::vec, arma::vec) \
                , arma::mat *X \
                , int N \
                , int M \
                , double dT0 \
                , double dStep) 
{
    // Function: vFuncIn( t, X(i), X(i+1) )

    // Routine
    double ti = dT0;
    
    for (int iii=0; iii<N-1; iii++) {
        // Iterate until the end of the array
        for (int j=0; j<M; j++) {
            // Iterate through each equation
            // Recalculate each of them as dStep could change (later) 
            X->col(iii+1) = X->col(iii) + \
                          (*vFuncIn)(ti, X->col(iii), X->col(iii+1)) *dStep;
        }
        ti += dStep;

    }


}


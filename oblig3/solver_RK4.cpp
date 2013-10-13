// RUNGE-KUTTA 4th order SOLVER ROUTINE

// Swallows functions, initial values and solution arrays

//#include <iostream>
//#include <cmath>
#include <armadillo>

#include "headers.h"

using namespace std;
using namespace arma;

void solver_RK4(   arma::vec (*vFuncIn)(arma::vec, arma::vec) \
                , arma::mat *X \
                , int N \
                , int M \
                , double dStep) 
{
    // Function: vFuncIn( X(i), X(i+1) )

    // Routine
    
    for (int iii=0; iii<N-1; iii++) {
        // Iterate until the end of the array
        for (int j=0; j<M; j++) {
            // Iterate through each equation
            // Recalculate each of them as dStep could change (later)
            
            arma::vec k1 = (*vFuncIn)(X->col(iii), X->col(iii+1));
            arma::vec k2 = (*vFuncIn)(X->col(iii) + k1*0.5*dStep,\
                   X->col(iii) + k1*0.5*dStep) ;
            arma::vec k3 = (*vFuncIn)(X->col(iii) + k2*0.5*dStep,\
                   X->col(iii) + k2*0.5*dStep) ;
            arma::vec k4 = (*vFuncIn)(X->col(iii) + k3*dStep, \
                   X->col(iii) + k3*dStep) ;



            X->col(iii+1) = X->col(iii) + \
                          1./6 * dStep * (k1 + 2*k2 + 2*k3 + k4);
        }
    }


}


// RUNGE-KUTTA 4th order SOLVER ROUTINE

// Swallows functions, initial values and solution arrays

#include <iostream>
//#include <cmath>
#include <armadillo>

#include "headers.h"

using namespace std;
using namespace arma;

void solver_RK4(   arma::vec (*vFuncIn)(double, arma::vec) \
                , arma::mat *X \
                , int N \
                , int M \
                , double dT0 \
                , double dStep) 
{
    // Function: vFuncIn( t, X(i) )

    // Routine
    
    arma::vec k1(M), k2(M), k3(M), k4(M);
    double ti = dT0;

    for (int iii=0; iii<N-1; iii++) {
        // Iterate until the end of the array
            
        k1 = (*vFuncIn)(ti, X->col(iii));
        k2 = (*vFuncIn)(ti + 0.5*dStep,\
               X->col(iii) + k1*0.5*dStep) ;
        k3 = (*vFuncIn)(ti + 0.5*dStep,\
               X->col(iii) + k2*0.5*dStep) ;
        k4 = (*vFuncIn)(ti + dStep, \
               X->col(iii) + k3*dStep) ;

        X->col(iii+1) = X->col(iii) + \
                      1./6 * dStep * (k1 + 2*k2 + 2*k3 + k4);

        ti += dStep;
    }


}


// EULER CROMER SOLVER ROUTINE

// Swallows functions, initial values and solution arrays

#include <iostream>
#include <cmath>
#include <armadillo>

#include "headers.h"

using namespace std;
using namespace arma;

void solver_EC(   arma::vec (*vFuncIn)(arma::vec) \
                , arma::vec *X0 \
                , arma::mat *X \
                , int N) 
{
    vec X2(3);
        X2 = (*vFuncIn)(*X0);

    X2.print();

}


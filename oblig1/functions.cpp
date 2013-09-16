// SOLVER FUNCTION
//
// Method: 
//  1) iterating forward, applying row multiplication factors (pivoting)
//  2) iterating backward, applying row m. factors (back substitution)
//  3) iterating forward, applying normalisation factor (gives leading 1s)
#include <iostream>
#include "time.h"

void solver(double *fX, double fH, double N)
{
    //begin timer
    clock_t start2, finish2;
    start2 = clock();
    
    
    // 1) Forward: 
    for (int iii=1; iii<N; iii++) {
        fX[iii] += iii/(iii + 1.0)  * fX[iii-1];
    }
    
    // 2) Backward:
    for (int jjj=N-2; jjj>= 0; jjj--) {
        fX[jjj] += (jjj + 1.0)/(jjj+2.0) * fX[jjj+1];
    }

    // 3) Normalise:
    for (int kkk=0; kkk<N; kkk++) {
        fX[kkk] *= (kkk+1.0)/(kkk+2.0);
    }


    // stop timer
    finish2 = clock();

    std::cout << "Time elapsed: tridiagonal solution scheme:" << std::endl;
    std::cout << ((finish2 - start2) / CLOCKS_PER_SEC ) << std::endl;
}

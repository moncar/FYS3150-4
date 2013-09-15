// SOLVER FUNCTION
//
// Method: 
//  1) iterating forward, applying row multiplication factors (pivoting)
//  2) iterating backward, applying row m. factors (back substitution)
//  3) iterating forward, applying normalisation factor (gives leading 1s)

#include <iostream>
#include <math.h>
using namespace std;
void solver(double *fX, double fH, double N)
{
    double j = 0;
    
    // 1) Forward: 
    for (int iii=1; iii<N; iii++) {
        j = iii + 1.0; 
        fX[iii] += (j-1)/j * fX[iii-1];
    }
    
    // 2) Backward:
    for (int jjj=N-2; jjj>= 0; jjj--) {
        j = jjj+1;
        fX[jjj] += (j+1)/(j+2) * fX[jjj+1];
    }

    // 3) Normalise:
    for (int kkk=0; kkk<N; kkk++) {
        j = kkk+1;
        fX[kkk] = pow(fH,2) *  j/(j+2) * fX[kkk];
    }

}

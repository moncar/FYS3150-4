// Jacobi eigenvalue algorithm
// Rotations until
//  - offdiagonal elements are zero,
//  - diagonal elements are eigenvalues
//

#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

void jacobiRotation(arma::mat &A, int N) {
    // Takes a sym. NxN matrix A and rotates it
    
    double eps  = 1e-3;     // precision
    int n       = 0;        // counter
    double offDiag=eps+1;

    while (offDiag > eps && n < N*N) {
        cout << n << endl;

        // Find largest off-diagonal element -> needed to find angle theta
        mat Aoff =  A; //zeros<mat>(N,N);
        Aoff.diag() = zeros<vec>(N);
        //Aoff.diag(+1) = A.diag(+1);
        //Aoff.diag(-1) = A.diag(-1);   // omitting main diagonal

        uword rowN; uword colN;
        int k; int l;
        double maxA = Aoff.max(rowN,colN);

        // Find sum over all offdiagonal elements,
        // % operator: elementwise multiplication
        // accu      : accumulates all elements
        offDiag = accu(Aoff % Aoff);
        cout << "Offdiag = " << offDiag << endl;

        if (maxA <= 0) {
            // first run -> chose random index,
            k = 0; l = 1;
            cout << A(k,l) << " at " <<  k << l << endl;
        }
        else {
            k = rowN > colN ? rowN : colN;
            l = rowN < colN ? rowN : colN;
            cout << maxA << " at " <<  rowN << colN << endl;
        }


        // find angle theta

        double cot2Theta = (A(l,l) - A(k,k))/(2*A(k,l));
        double tanTheta  = min((-cot2Theta + sqrt(1 + pow(cot2Theta,2))), \
                (-cot2Theta + sqrt(1 + pow(cot2Theta,2))) );

        // find: cos(theta) and sin(theta) from tan(theta)

        double cosTheta  = 1./ (sqrt( 1 + pow(tanTheta,2)));
        double sinTheta  = tanTheta * cosTheta;

        
        // Apply rotation
        
       for (int j = 0; j<N; j++) {
           if (j != k && j != l) {
               A(j,k) = A(j,k) * cosTheta - A(j,l) * sinTheta;
               A(j,l) = A(j,l) * cosTheta + A(j,k) * sinTheta;
           }}
           A(k,k) = A(k,k)*cosTheta*cosTheta \
                    - 2*A(k,l) * cosTheta * sinTheta \
                    +   A(l,l) * sinTheta * sinTheta;
           A(l,l) = A(l,l) * cosTheta*cosTheta \
                    + 2*A(k,l) * cosTheta * sinTheta \
                    +   A(k,k) * sinTheta * sinTheta;
           A(k,l) = (A(k,k) - A(l,l)) * cosTheta * sinTheta \
                    +   A(k,l) * (cosTheta * cosTheta - sinTheta *sinTheta);
        n += 1;
    }
}

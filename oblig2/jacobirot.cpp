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
    
    // Find largest off-diagonal element -> needed to find angle theta
    mat Aoff =  A; //zeros<mat>(N,N);
    //Aoff.diag(+1) = A.diag(+1);
    //Aoff.diag(-1) = A.diag(-1);   // omitting main diagonal

    while (offDiag > eps && n < N*N) {
        cout << n << endl;

        Aoff = A;
        Aoff.diag() = zeros<vec>(N);

        uword rowN; uword colN;
        int k; int l;
        double maxA = Aoff.max(rowN,colN);

        // Find sum over all offdiagonal elements,
        // % operator: elementwise multiplication
        // accu      : accumulates all elements
        offDiag = sqrt(accu(Aoff % Aoff));
        cout << "Offdiag = " << offDiag << endl;

        if (maxA <= 0) {
            // first run -> chose random index,
            k = 0; l = 1;
            cout << A(k,l) << " at " <<  k << l << endl;
        }
        else {
            // find 1 < k < l < n
            k = rowN < colN ? rowN : colN;
            l = rowN > colN ? rowN : colN;
            cout << k << l << endl;
            cout << maxA << " at A(" << rowN+1<< ","<< colN+1<< ")" << endl;
        }


        // find angle theta

        double cot2Theta = (A(l,l) - A(k,k))/(2*A(k,l));
        
        double tanTheta  = min(( -1./(-cot2Theta - sqrt(1 + pow(cot2Theta,2)))), ( -1./(-cot2Theta + sqrt(1 + pow(cot2Theta,2)))) ) ;

        cout << tanTheta << endl;
        cout << cot2Theta << "\t" << sqrt(1+ pow(cot2Theta,2)) << endl;
        
        // find: cos(theta) and sin(theta) from tan(theta)

        double cosTheta  = 1./ (sqrt( 1 + pow(tanTheta,2)));
        double sinTheta  = tanTheta * cosTheta;

        
        // Apply rotation
        
       double Akk = A(k,k);
       double All = A(l,l);
       double Akl = A(k,l);
       for (int j = 0; j<N; j++) {
           double Ajk = A(j,k);
           double Ajl = A(j,l);
           if (j != k && j != l) {
               A(j,k) = Ajk * cosTheta - Ajl * sinTheta;
               A(j,l) = Ajl * cosTheta + Ajk * sinTheta;
           }}

       A(k,k) = Akk*cosTheta*cosTheta \
                - 2*Akl * cosTheta * sinTheta \
                +   All * sinTheta * sinTheta;
       A(l,l) = All * cosTheta*cosTheta \
                + 2*Akl * cosTheta * sinTheta \
                +   Akk * sinTheta * sinTheta;
       A(k,l) = (Akk - All) * cosTheta * sinTheta \
                +   Akl * (cosTheta * cosTheta - sinTheta *sinTheta);
       
        n += 1;
    }
}

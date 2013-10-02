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
        
    double eps  = 1e-8;     // precision
    int n       = 0;        // counter
    double offDiag=eps+1;
    
    // Matrix to hold off-diagonal values 
    mat Aoff =  A; //zeros<mat>(N,N);
    //Aoff.diag(+1) = A.diag(+1);
    //Aoff.diag(-1) = A.diag(-1);   // omitting main diagonal
     
    // Matrix to hold absolute values of off-diagonal elements
    mat AoffAbs = A;

    while (offDiag > eps && n < N*N*N) {
        // Break if:
        // - The sum of offdiagonal elements are smaller than a epsilon
        // - iteration no n is above a limit
        cout << n << "========= iteration =========" << n << endl;

        Aoff = A;
        Aoff.diag() = zeros<vec>(N);

        // Find maximum of absolute values of off-diagonal elements:
        uword rowN; uword colN;
        int k; int l;
        AoffAbs = abs(Aoff);
        double maxA = AoffAbs.max(rowN, colN);

        //A.print();
        cout << "\n" << endl;
        //Aoff.print();

        // Find sum over all offdiagonal elements,
        // % operator: elementwise multiplication
        // accu      : accumulates all elements
        offDiag = sqrt(accu(Aoff % Aoff));
        cout << "Offdiag = " << offDiag << endl;

        // find 1 < k < l < n <-- ERROR in compendium?
        k = rowN; //rowN < colN ? rowN : colN;
        l = colN; //rowN > colN ? rowN : colN;
        //cout << "k=" <<  k << " l=" << l << endl;
        cout << maxA << " at A(" << rowN+1<< ","<< colN+1<< ")" << endl;


        // find angle theta
        //cout << "A(l,l)=" << A(l,l) << " A(k,k)=" << A(k,k) << " A(k,l)=" << A(k,l) << endl;
        //cout << "Finding angle theta: printing tanTheta and cot2Theta" << endl;
          //
        double cot2Theta = (A(l,l) - A(k,k))/(2*A(k,l));

        double tanTheta  = min( (-cot2Theta + sqrt(1+ cot2Theta*cot2Theta)) \
                                            -cot2Theta - sqrt(1+ cot2Theta*cot2Theta) );
              
              //= min(( -1./(-cot2Theta - sqrt(1 + pow(cot2Theta,2)))), ( -1./(-cot2Theta + sqrt(1 + pow(cot2Theta,2)))) ) ;

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
                //cout << j<<k<<l << endl;
                A(j,k) = Ajk * cosTheta - Ajl * sinTheta;
                A(j,l) = Ajl * cosTheta + Ajk * sinTheta;
                //cout << "A(j,k)=" << A(j,k) << endl;
                //cout << "A(j,l)=" << A(j,l) << endl;
           }
        }

        A(k,k) = Akk*cosTheta*cosTheta \
                - 2*Akl * cosTheta * sinTheta \
                +   All * sinTheta * sinTheta;
        //cout << "A(k,k)=" << A(k,k) << endl;
        A(l,l) = All * cosTheta*cosTheta \
                + 2*Akl * cosTheta * sinTheta \
                +   Akk * sinTheta * sinTheta;
        //cout <<  "A(l,l)=" <<A(l,l) << endl;
        A(k,l) = (Akk - All) * cosTheta * sinTheta \
                +   Akl * (cosTheta * cosTheta - sinTheta *sinTheta);
        //cout <<  "A(k,l)=" <<A(k,l) << endl;

        n += 1;
    }
}

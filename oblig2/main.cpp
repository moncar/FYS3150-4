// 2013-10-07 FYS4150
// Marius Berge Eide

// Project 2
// Solving the Schrodinger equation

#include <iostream>
#include <cmath>
#include <armadillo>

#include "headers.h"
#include "lib.h"

using namespace std;
using namespace arma;

// Constants

const double hbar = 6.58211928e-16;   // eV s
const double m_e  = 0.51100e6;        // eV c-2

inline double fPot1(double fRho_min, double fRho_h, int i) {
    return pow((fRho_min + i*fRho_h), 2);
}

inline double fPot2(double fRho_min, double fRho_h, int i, double omega) {
    return pow(( omega * (fRho_min + i*fRho_h) ), 2) \
            + 1./(fRho_min + max(i,1)*fRho_h);
    // using max(i,1) to prevent division w/0
}

int main(int argc, char* argv[]) {

    // Physics: potential of the harmonic oscillator V(r)
    // Potential: V(rho) = rho^2
    int N;
    if (argc > 1)   { N = atoi(argv[1]); }
    else            { N = 12; }
    double fRho_min = 0.0;
    double fRho_max = 100.0; // 1./fAlpha * r // 10 alphas
    double omega    = 0.01; // oscillator "frequency"

    // Arrays to hold eigenvalues
    vec eigvals1(N);
    vec eigvals2(N);
    vec eigvals3(N);

    cout << "JACOBI" << endl;
    vStartJacobi(N, fRho_min, fRho_max, omega, &eigvals1);
    cout << "ARMADILLO" << endl;
    vStartEIGVAL(N, fRho_min, fRho_max, omega, &eigvals2);
    cout << "TQLI HOUSEHOLDER" << endl;
    vStartTQLI  (N, fRho_min, fRho_max, omega, &eigvals3);

    // Pass on 3*N long array with results to file, to be stored as cols
    // Needs: rho values

    double* rhos = new double[3*N];
    double* vals = new double[3*N];

    for (int j=0; j < 3; j++) {
        // three columns as one
        for (int i=0; i < N; i++) {
            rhos[i + j*N ] = fRho_min + i*(fRho_max-fRho_min)/N;
            if (j==0) { vals[i]         = eigvals1(i); } 
            if (j==1) { vals[i + j*N]   = eigvals2(i); } 
            if (j==2) { vals[i + j*N]   = eigvals3(i); }
        }
    }

    if (argc <= 2) {
        cout << "ERROR! Missing output file specification:\n \
                 \t obl2.x N outputfile.txt" << endl;
        exit(1);
    }
    else {
        outputFile(N, 3, rhos, vals, &argv[2]); 
    }

    return 0;

}

void vStartJacobi( \
        int N, \
        double fRho_min, \
        double fRho_max, \
        double omega, \
        arma::vec* eigvals){
    // Function to initialise and run Jacobi rotation
    // Creates tridiagonal matrix and passes it on to jacobirot.cpp

    double fRho_h   = (fRho_max - fRho_min)/N;
    //double fOmega   = 
    //double fK       = m_e * pow(fOmega, 2);
    //double fAlpha   = pow((pow(hbar, 2)/(m_e * fK)), 1./4);
    
    // Set up a potential, V and diagonal elemnts d
    vec fV(N);  // potential
    vec fD(N);  // diagonal
    
    for (int i=0; i<N; i++) {
        // V_i   = rho_i ^ 2
        // rho_i = rho_min + i*h
        fV[i] = fPot1(fRho_min, fRho_h, i);     // potential
        //fV[i] = fPot2(fRho_min, fRho_h, i, omega); // 2 electrons
        fD[i] = fV[i] + 2.*pow(fRho_h, -2);     // diagonal
    }

    // Set up a tridiagonal matrix
    mat A       = zeros<mat>(N,N);

    A.diag()    = fD;
    A.diag(+1)  -= pow(fRho_h, -2.); 
    A.diag(-1)  -= pow(fRho_h, -2.); 

    A.print();

    // Apply operations on A to find eigenvalues
    jacobiRotation(A, N);

    A.print();

    // Eigenvalues are on the diagonal
    *eigvals = diagvec(A);
    *eigvals = sort(*eigvals);

    //cout << &A << endl;
    //cout << sort(*eigvals) << endl;

}


void vStartEIGVAL(\
        int N, \
        double fRho_min, \
        double fRho_max, \
        double omega, \
        arma::vec* eigvals ){
    // Function to initialise tridiagonal arrays
    // and pass these on to EIG_SYM method of Armadillo
    // Prints eigenvalues and eigenvectors

    double fRho_h   = (fRho_max - fRho_min)/N;
    //double fOmega   = 
    //double fK       = m_e * pow(fOmega, 2);
    //double fAlpha   = pow((pow(hbar, 2)/(m_e * fK)), 1./4);
    
    // Set up a potential, V and diagonal elemnts d
    vec fV(N);  // potential
    vec fD(N);  // diagonal
    
    for (int i=0; i<N; i++) {
        // V_i   = rho_i ^ 2
        // rho_i = rho_min + i*h
        fV[i] = fPot1(fRho_min, fRho_h, i);     // potential
        //fV[i] = fPot2(fRho_min, fRho_h, i, omega); // 2 electrons
        fD[i] = fV[i] + 2.*pow(fRho_h, -2);     // diagonal
    }

    // Set up a tridiagonal matrix
    mat A       = zeros<mat>(N,N);

    A.diag()    = fD;
    A.diag(+1)  -= pow(fRho_h, -2.); 
    A.diag(-1)  -= pow(fRho_h, -2.); 

    A.print();

    // Apply operations on A to find eigenvalues
    //vec eigvals;
    mat eigvec;

    eig_sym(*eigvals, eigvec, A);

    eigvals->print();
    // Woho! Same as 
    // (*eigvals).print();
}


void vStartTQLI(\
        int N, \
        double fRho_min, \
        double fRho_max, \
        double omega, \
        arma::vec* eigvals ) {
    // Calling TQLI: Householder's method for finding eigenvalues and
    // eigenvectors.

    double fRho_h   = (fRho_max - fRho_min)/N;
    //double fOmega   = 
    //double fK       = m_e * pow(fOmega, 2);
    //double fAlpha   = pow((pow(hbar, 2)/(m_e * fK)), 1./4);
    
    // Set up a potential, V and diagonal elemnts d
    double* fV2= new double[N];   // potential
    double* fD2= new double[N];   // diagonal
    double* fO2= new double[N];   // off-diagonal elements
    // C++ matrix to hold eigenvectors
    double** fEV = new double*[N];
    for (int i=0; i<N; i++) {
        fEV[i] = new double[N];
    }
    
    for (int i=0; i<N; i++) {
        // V_i   = rho_i ^ 2
        // rho_i = rho_min + i*h
        fV2[i] = fPot1(fRho_min, fRho_h, i);     // potential
        //fV2[i] = fPot2(fRho_min, fRho_h, i, omega); // 2 electrons
        fD2[i] = fV2[i] + 2.*pow(fRho_h, -2);       // diagonal
        //if (i<N-1) {
            fO2[i] = pow(fRho_h, -2.); 
        //}
    }

    // Apply Householder's algorithm to find eigenvalues
    tqli(fD2, fO2, N, fEV);

    // Sort eigenvalues
    sort(fD2, fD2+N);

    for (int i=0; i<N; i++) {
        (*eigvals)(i) = fD2[i];
        cout << fD2[i] << endl;
    }

    // Housekeeping
    delete [] fV2;
    delete [] fD2;
    delete [] fO2;
    delete [] fEV;
}




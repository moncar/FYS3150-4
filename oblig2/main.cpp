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

int main(int argc, char* argv[]) {

    // Physics: potential of the harmonic oscillator V(r)
    // Potential: V(rho) = rho^2
    int N = 10;
    double fRho_min = 0.0;
    double fRho_max = 10.0; // 1./fAlpha * r // 10 alphas

    vStartJacobi(N, fRho_min, fRho_max);
    vStartEIGVAL(N, fRho_min, fRho_max);
    vStartTQLI(N, fRho_min, fRho_max);

    return 0;

}

void vStartJacobi(double N, double fRho_min, double fRho_max) {
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
        fV[i] = pow((fRho_min + i*fRho_h), 2);  // potential
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
    vec eig = diagvec(A);

    //cout << &A << endl;
    cout << sort(eig) << endl;

}


void vStartEIGVAL(double N, double fRho_min, double fRho_max) {
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
        fV[i] = pow((fRho_min + i*fRho_h), 2);  // potential
        fD[i] = fV[i] + 2.*pow(fRho_h, -2);     // diagonal
    }

    // Set up a tridiagonal matrix
    mat A       = zeros<mat>(N,N);

    A.diag()    = fD;
    A.diag(+1)  -= pow(fRho_h, -2.); 
    A.diag(-1)  -= pow(fRho_h, -2.); 

    A.print();

    // Apply operations on A to find eigenvalues
    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, A);

    eigval.print();
}


void vStartTQLI(int N, double fRho_min, double fRho_max) {
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
        fV2[i] = pow((fRho_min + i*fRho_h), 2);     // potential
        fD2[i] = fV2[i] + 2.*pow(fRho_h, -2);       // diagonal
        //if (i<N-1) {
            fO2[i] = pow(fRho_h, -2.); 
        //}
    }

    // Apply Householder's algorithm to find eigenvalues
    tqli(fD2, fO2, N, fEV);

    for (int i=0; i<N; i++) {
        cout << fD2[i] << endl;
    }

    // Housekeeping
    delete [] fV2;
    delete [] fD2;
    delete [] fO2;
    delete [] fEV;
}




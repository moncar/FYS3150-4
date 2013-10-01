// 2013-10-07 FYS4150
// Marius Berge Eide

// Project 2
// Solving the Schrodinger equation

#include <iostream>
#include <cmath>
#include <armadillo>

#include "headers.h"

using namespace std;
using namespace arma;

// Constants

const double hbar = 6.58211928e-16;   // eV s
const double m_e  = 0.51100e6;        // eV c-2

int main(int argc, char* argv[]) {

    // Physics: potential of the harmonic oscillator V(r)
    // Potential: V(rho) = rho^2
    int N = 25;
    double fRho_min = 0.0;
    double fRho_max = 10.0; // 1./fAlpha * r    // 10 alphas
    double fRho_h   = (fRho_max - fRho_min)/N;
    //double fOmega   = 
    //double fK       = m_e * pow(fOmega, 2);
    //double fAlpha   = pow((pow(hbar, 2)/(m_e * fK)), 1./4);
    
    // Set up a potential, V and diagonal elemnts d
    vec fV(N);
    vec fD(N);
    
    for (int i=0; i<N; i++) {
        // V_i   = rho_i ^ 2
        // rho_i = rho_min + i*h
        fV[i] = pow((fRho_min + i*fRho_h), 2);
        fD[i] = fV[i] + 2.*pow(fRho_h, -2);
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

    return 0;

}

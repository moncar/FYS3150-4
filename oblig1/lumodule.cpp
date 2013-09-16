// LU part of project 1
//

#include <armadillo>
#include <iostream>
#include "lib.h"
#include "time.h"

//using namespace arma;
using namespace std;


void LUmethod(int N, double* fSrcX)
// Input:   No. of iterations/size of matrix: N
//          N-sized vector holding source function
// Returns: Void. Changes fSrcX into approximation to
//          solution x in "Ax = b" where b is fSrcX.
{
    // Initiate tridiagonal matrix

    /*
    // Method A: using Armadillo
    mat A = zeros<mat>(N,N);
    A.diag()   += 2;
    A.diag(-1) -= 1;
    A.diag(1)  -= 1;

    std::cout << A[0] << std::endl;

    // Method B: using hardcore C++
    double **B = new double*[N];
    for (int i=0; i<N; i++) {
        B[i] = new double[N];
    }

    std::cout << "MATRIX B: " << N <<  std::endl;
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if (i == j) { B[i][j] = 2.0; }
            else if (i != j && abs(i-j) < 2) { B[i][j] = -1.0; }
            else {B[i][j] = 0.0;}
            std::cout << B[i][j] << "\t";
            if ((j+1) % 5 == 0) { std::cout << "\n"; }
        }
    }

    */


    // Method C: using Blitz implementation in lib.h
    double **C;
    C = (double **) matrix(N,N, sizeof(double));
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            if (i == j) { C[i][j] = 2.0; }
            else if (i != j && abs(i-j) < 2) { C[i][j] = -1.0; }
            else {C[i][j] = 0.0;}
            //std::cout << C[i][j] << "\t";
            //if ((j+1) % N == 0) { std::cout << "\n"; }
        }
    }


    // LU decompose matrix
    
    //begin timing
    clock_t start, finish;
    start = clock();

    
    // ludcmp(double a, int n, int indx, double d)
    // where    a: NxN matrix -> LU-factorised -> becomes LU
    //          n: dimensionality
    //       indx: index vector, holds row permutations
    //          d: outputs +1 or -1: number of interchanges are even or odd
    int *indx = new int[N];
    double d = 0.0;

    ludcmp(C, N, indx, &d);

    // Print matrix
    /*
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            std::cout << C[i][j] << "\t";
            if ((j+1) % N == 0) { std::cout << "\n"; }
        }
    } */

    // Solve for x in "Ax = b", using input array fSrcX as b

    lubksb(C, N, indx, fSrcX);

    // stop timer
    finish = clock();

    std::cout << "Time elapsed: LU solution scheme:" << std::endl;
    cout << ((finish - start) / CLOCKS_PER_SEC ) << endl;
} 

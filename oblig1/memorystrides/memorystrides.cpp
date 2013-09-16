#include <iostream>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
    // Multiplication of two matrices, A = BC. 
    
    int N = 1000;

    double **A = new double*[N];
    double **B = new double*[N];
    double **C = new double*[N];

    for (int i=0; i < N; i++) {
        A[i] = new double[N];
        B[i] = new double[N];
        C[i] = new double[N];
    }


    // Initialise timers
    clock_t start1, finish1;
    clock_t start2, finish2;
    clock_t start3, finish3;

    // Method 1: row major order iterations
    start1 = clock();
    for (int i=0; i < N; i++) {
        for (int j=0; j < N; j++) {
            for (int k=0; k < N; k++) {
                A[i][j] += B[i][k] * C[k][j];
            }
        }
    }
    finish1 = clock();
    cout << "Time elapsed: Row major order iterations in matrix multiplication: (in seconds):" << endl;
    cout << ((finish1 - start1) / CLOCKS_PER_SEC ) << endl;

    // Method 2: column major order iterations
    start2 = clock();
    for (int j=0; j< N; j++) {
        for (int i=0; i < N; i++) {
            for (int k=0; k < N; k++) {
                A[i][j] += B[i][k] * C[k][j];
            }
        }
    }
    finish2 = clock();
    cout << "Time elapsed: Column major order iterations in matrix multiplication: (in seconds):" << endl;
    cout << ((finish2 - start2) / CLOCKS_PER_SEC ) << endl;


    // Method 3: Armadillo/BLAS multiplication

    start3 = clock();

    mat Am = zeros<mat>(N,N);
    mat Bm = randu<mat>(N,N);
    mat Cm = randu<mat>(N,N);

    Am = Bm * Cm;

    finish3 = clock();

    cout << "Time elapsed: Armadillo matrix multiplication: (in seconds):" << endl;
    cout << ((finish3 - start3) / CLOCKS_PER_SEC ) << endl;


    delete [] A;
    delete [] B;
    delete [] C;


    return 0;
}

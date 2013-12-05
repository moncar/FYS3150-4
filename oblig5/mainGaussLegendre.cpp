#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lib.h"
#include "headers.h"

int mainGaussLegendre(int argc, char* argv[], int N) {
    
    cout << "Launching Gauss-Legendre integrator. N=" << (int)N << endl;

    // Gauss-Legendre:  integrate wave-function
    //                  which is approx zero at r=5

    double lima = -5.0;
    double limb = +5.0;


    double* x = new double[N];  // zeros/ mesh points
    double* w = new double[N];  // weights
    double intsum = 0.0;        // integration sum

    // init: mesh points and weights
    gauleg(lima, limb, x,w, N);

    // integrate
    #pragma omp parallel for reduction(+:intsum)
    for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
    for (int k=0; k < N; k++) {
    for (int l=0; l < N; l++) {
    for (int m=0; m < N; m++) {
    for (int n=0; n < N; n++) {
        intsum += w[i]*w[j]*w[k]*w[l]*w[m]*w[n] \
                * wavefunc(x[i], x[j], x[k], x[l], x[m], x[n]);
    }}}}}}

    printf("Integral sum: %g , with N=%d \n", intsum, N);
    
    // Housekeeping
    delete [] x;
    delete [] w;
    
    
    return 0;
}

#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lib.h"
#include "headers.h"

int mainGaussHermite(int argc, char* argv[], int N) {
    
    cout << "Launching Gauss-Hermite integrator. N=" << (int)N << endl;

    // Gauss-Hermite:   integrate wave-function
    //                  using ortho func: Hermite, +/- infty
    
    double* xx = new double[N+1]; // zeros/mesh points
    double* ww = new double[N+1]; // weights
    double intsum2= 0.0;        // integration sum: G-H

    // init: mesh points and weights
    GaussHermite(xx, ww, N);

    cout << ww[0] << " " << ww[1] << endl;

    // integrate
    #pragma omp parallel for reduction(+:intsum2)
    for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
    for (int k=0; k < N; k++) {
    for (int l=0; l < N; l++) {
    for (int m=0; m < N; m++) {
    for (int n=0; n < N; n++) {
        intsum2 += ww[i]*ww[j]*ww[k]*ww[l]*ww[m]*ww[n] \
                * wavefunc_red(xx[i], xx[j], xx[k], xx[l], xx[m], xx[n]);
    }}}}}}

    printf("Integral sum: %.12g , with N=%d \n", intsum2, N);
    
    // Housekeeping
    delete [] xx;  
    delete [] ww;

    return 0;

}


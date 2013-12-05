#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lib.h"
#include "headers.h"

int mainMonteCarloBruteForce(int argc, char* argv[], int N) {
    
    // Monte-Carlo
    
    cout << "Launching Monte Carlo integrator. N=" << N << endl;


    // Brute-force Monte Carlo

    // init random seed 
    //srand(time(0));
    
    // Limits:
    double limit = 5.0;

    // GNU Scientific library implementation
    const gsl_rng_type * T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;    // generator function
    r = gsl_rng_alloc(T);   // random number generator

    // C++2011 implementation:
    // random_device rd;
    // uniform_real_distribution<double> uniform_real(-5.0,5.0);
   
    // Jacobi determinant
    // length[-a,a] = 2a
    double jacobiDet  = pow((2*limit), 6.0);

    double** xxx = new double*[6];
    double*  fx  = new double[N];
    #pragma omp parallel for
    for (int i=0; i < 6; i++) {
        xxx[i] = new double[N];
    }

    for (int j=0; j < 6; j++) {
        #pragma omp parallel for
        for (int i=0; i < N; i++) {
            xxx[j][i]   =  gsl_ran_flat(r, -limit, +limit); 
                        //  C++2011:
                        //= uniform_real(rd);
        }
    }


    double intsum3 = 0.0;
    double sigmasq = 0.0;
    double variance= 0.0;

    for (int j=0; j < 6; j++) {
        cout << xxx[j][0] << "  \t" << xxx[j][1] << endl;
    }

    // Calculate f(x) at the given points
    #pragma omp parallel for
    for (int i=0; i < N; i++) {
        fx[i]    = wavefunc(xxx[0][i], xxx[1][i], xxx[2][i] \
                            , xxx[3][i], xxx[4][i], xxx[5][i]);
    }

    #pragma omp parallel for reduction(+:intsum3,sigmasq)
    for (int i=0; i < N; i++) {
        intsum3 += fx[i];
        sigmasq += fx[i]*fx[i];
    }

    double integral = intsum3/N;
    sigmasq         = sigmasq/N;
    variance        = sigmasq - integral*integral;
    variance        = variance / N ; // ???

    // Scale
    integral *= jacobiDet;
    double stddev = sqrt(variance) * jacobiDet;

    printf("Integral sum: %.12g , with N=%d \n", integral, N);
    cout << "\t" << sigmasq/N << " x^2: " << integral*integral/jacobiDet << endl;
    printf("\tVariance: %.12g,\n\tstd.dev: %.12g\n", variance, stddev);

    if (argc > 3) { 
        outputFile2(N, 1, xxx[1], fx, &argv[3]);
    }
    
    // Housekeeping
    gsl_rng_free (r);
    delete [] xxx;
    delete [] fx;

    return 0;
}

#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lib.h"
#include "headers.h"

int mainMonteCarloImportanceSampling(int argc, char* argv[], int N) {

    // Monte Carlo: Importance sampling

    // Need: - timer for all methods
    //       - variance and stddev for all methods

    // Initialise
    double** x4 = new double*[6];   // holding: random numbers
    double*  fx4= new double[N];    // holding: function values
    clock_t t;                      // timer
    t = clock();

    // Set std.dev for Gaussian function
    double sigma = 1./sqrt(2);
    // Set Jacobi-det for integral, following MHJ
    double jacobiDet = pow(4*atan(1), 3);


    // Create random numbers and store in array x4 w/ 6 dimensions
    RandNoGenGauss(x4, sigma, 6, N);

    double intsum4 = 0.0;
    double sigmasq = 0.0;
    double variance= 0.0;

    for (int j=0; j < 6; j++) {
        cout << x4[j][0] << "  \t" << x4[j][1] << endl;
    }

    // Calculate f(x) at the given points
    #pragma omp parallel for
    for (int i=0; i < N; i++) {
        fx4[i]   = wavefunc_red(x4[0][i], x4[1][i], x4[2][i] \
                            , x4[3][i], x4[4][i], x4[5][i]);
    }

    cout << "Done calculating wavefunc" << endl;

    #pragma omp parallel for reduction(+:intsum4,sigmasq)
    for (int i=0; i < N; i++) {
        intsum4 += fx4[i];
        sigmasq += fx4[i]*fx4[i];
    }

    cout << "Done w/ reduction" << endl;

    double integral = intsum4/N;
    sigmasq         = sigmasq/N;
    variance        = sigmasq - integral*integral;
    variance        = variance / N ; // ???

    // Scaling
    integral        = integral * jacobiDet;
    double stddev   = sqrt(variance) * jacobiDet;

    printf("Integral sum: %.12g , with N=%d \n", integral, N);
    printf("\tStd.dev: %.12g\n", stddev);
    

    // Time spent?
    t = clock() - t;
    cout << "Elapsed time: " << (float)t/CLOCKS_PER_SEC << "s" << endl;

    // Save if argument says so
    if (argc > 3) { 
        outputFile2(N, 1, x4[0], fx4, &argv[3]);
    }

    // Housekeeping

    delete [] x4;
    delete [] fx4;

    return 0;    

}

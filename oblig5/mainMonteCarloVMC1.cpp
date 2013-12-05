#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lib.h"
#include "headers.h"

int mainMonteCarloVMC1(int argc, char* argv[], int N) {
    // Variational part: using test wave function PSI_T1
    // Solution strategy:
    //  - Create 6 random nos - one for each dim
    //  - Calc wave func
    //  - Check ratio: new wavefunc sq / old wavefunc sq >= ran[0,1]?
    //  - YES:  use random nos func as starting point,
    //          save wavefunc value,
    //          save random nos as locations,
    //          iteration++
    //  - NO: keep old ones, try new sample, iteration = iteration

    // Initialise
    double** xVMC1 = new double*[N];   // holding: random numbers
    double*  fxVMC1= new double[2*N];    // holding: function values
    double*  ergLoc= new double[N];    // holding: local energy values
    clock_t t;                      // timer
    t = clock();

    // Set std.dev for Gaussian function
    double sigma = 1./sqrt(2);
    // Set Jacobi-det for integral, following MHJ
    double jacobiDet = pow(4*atan(1), 3);

    
    // Calculate random number
    // GNU Scientific library implementation
    const gsl_rng_type * T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_taus;    // generator function
    r = gsl_rng_alloc(T);   // random number generator
    
    // C++2011 implementation:
    //random_device generator;
    //normal_distribution<double> distribution(0.0, 2.0);

    // Initialise arrays
    #pragma omp parallel for
    for (int j=0; j < 6; j++) {
        xVMC1[j] = new double[N];
    }
    cout << "Done initialising." << endl;

    // Set initial trial wavefunc
    //

    //for (int j=0; j<6; j++) { xVMC1[j][0] = xi[j]; }


    
    int variations= 20;     // alpha variations
    double alphast= 0.01;    
    double alpha0 = 0.95;
    double* alphas= new double[variations];
    double* ergtot= new double[variations];
    double* stddev= new double[variations];
    double* variances = new double[variations];

    for (int a=0; a < variations; a++) { alphas[a] = alpha0 + a*alphast; }

    double  step = 0.0005;
    int     n;
    int     randcord;
    double  randno;
    double  randno2;
    double  randno3;
    double  probval;
    double  ratio;
    double  ergloc;
    double  ergloc_tmp;
    double  ergloc_sq;
    
    
    // Loop over all alphas
    #pragma omp parallel for private(ergloc,ergloc_tmp,ergloc_sq,n)
    for (int a=0; a< variations; a++) {
        cout << "Iterating for alpha=" << alphas[a] << endl;
        ergloc      = 0.0;
        ergloc_tmp  = 0.0;
        ergloc_sq   = 0.0;
        n       = 1;
        
        // Choose starting position
        double xi[6] = {1.0, 1.0, 1.0 \
                        ,2.0, 2.0, 2.0 };
        double xf[6] = {0.0, 0.0, 0.0 \
                    ,0.0, 0.0, 0.0 };
        fxVMC1[0] = wavefuncsqT1(xi, alphas[a]);


        while (n < N) {
            // Pick random coordinate [0,5] and calc ratio
            randcord = (int)floor(6*gsl_rng_uniform(r));
            // Calculate RVs
            randno = gsl_rng_uniform(r);    // in range [0,1).
            randno3= gsl_rng_uniform(r);    // in range [0,1).
            // Transform no -> [-1,1):
            randno2= 2*randno - 1.0;

            // Calculate new step:
            for (int j=0; j<6; j++) { xf[j] = xi[j]; }
            xi[randcord] = xf[randcord] + randno2*step;

            // Calc wavefunc ratio
            probval = wavefuncsqT1(xi, alphas[a]);
            ratio = probval \
                    /fxVMC1[n-1];

            // if ratio > randno
            // continue
            if (ratio > randno3) {
                for (int j=0; j<6; j++) { 
                    xVMC1[j][n] = xf[j];
                    xi[j] = xf[j];
                }
                // Calculate local energy
                ergloc_tmp= LocalEnergy1(xf[0], xf[1], xf[2] \
                                       , xf[3], xf[4], xf[5] \
                                       , alphas[a], N); 
                ergloc += ergloc_tmp;
                ergloc_sq += ergloc_tmp * ergloc_tmp;
                // Calculate probability value
                fxVMC1[n]= probval;
                n++;
            } 
        }
        // Calculate Hamiltonian 
        ergtot[a]   = 1./N * ergloc;

        // Calculate variance
        ergloc_sq   = 1./N * ergloc_sq;
        variances[a]= ergloc_sq - ergtot[a]*ergtot[a];
        stddev[a]   = sqrt(variances[a]/N);

        cout << "\tVar: " << variances[a] << ", ergloc_sq: " << ergloc_sq << endl;
        cout << "\tergloc: " << ergloc << endl;

        cout << "\tTotal energy: " << ergtot[a] << endl;
        cout << "\tStd.dev: " << stddev[a] << endl;
    }


    // Time spent?
    t = clock() - t;
    cout << "Elapsed time: " << (float)t/CLOCKS_PER_SEC << "s" << endl;
    
    double* xVMCvals = new double[N];
    for (int i=0; i<N; i++) { 
        xVMCvals[i] = xVMC1[0][i]; 
    }
//    for (int i=N; i<2*N; i++) {
//        fxVMC1[i] = ergLoc[i-N];
//    }

    if (argc > 3) { 
//        outputFile2(N, 1, xVMCvals, fxVMC1, &argv[3]);
        outputFile2(variations, 1, alphas, stddev, &argv[3]);        
    }

    // Housekeeping
    delete [] alphas;
    delete [] xVMC1;
    delete [] fxVMC1;
    delete [] xVMCvals;
    delete [] ergLoc;
    delete [] ergtot; 
    delete [] stddev;
    delete [] variances;

    return 0;    
}


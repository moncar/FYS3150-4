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
    // Timer
    double t = omp_get_wtime();     // get wall time

    //double (&wavefunc)(double*, double, double)  = wavefuncT1;
    //double (&wavefuncSq)(double*, double, double)= wavefuncT1sq;
    double (&wavefunc)(double*, double, double)  = wavefuncT2;
    double (&wavefuncSq)(double*, double, double)= wavefuncT2sq;

    // Variational MC:
    double alph0    =-1.2;
    double beta0    = 0.00;
    double alphStep = 0.005;
    double betaStep = 0.0078;
    int variations  = 128;

    // Variational matrix:
    //  - Alphas
    //  - Avg energies
    //  - Std.deviations
    double * alphas     = new double[variations];
    double * betas      = new double[variations];
    double * energies   = new double[variations];
    double **locergs    = new double*[variations];
    double * stddevs    = new double[variations];

    for (int a=0; a<variations; a++) {
        alphas[a] = alph0 + a*alphStep;
    }
    for (int b=0; b<variations; b++) {
        betas[b] = beta0 + b*betaStep;
    }
    for (int l=0; l<variations; l++) {
        locergs[l] = new double[2];
    }

    // Calculate random number
    // GNU Scientific library implementation
    const gsl_rng_type * T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_taus;       // generator function
    r = gsl_rng_alloc(T);   // random number generator

    // VARIATIONAL MONTE CARLO: VARIATE ALPHA & BETA
    cout << "Entering parallell section..." << endl;
    // In case: variate: ALPHA, fix BETA:
    int aplus=1;  int a = 0;
    int bplus=0;  int b = 0;
//betas[b] = 0.1;
    betas[b] = 0.086;
    
    // In case: variate: BETA, fix ALPHA:
//    int aplus=0;    int a = 0;
//    int bplus=1;    int b = 0;
//    alphas[a] = -1.0;

    // Chose smallest: number of avail threads or variations
    int numthr = min(omp_get_max_threads(), variations);
    cout << "No. parallel regions: " << numthr << endl;
    omp_set_num_threads(numthr);

    #pragma omp parallel for
    for (int varp=0; varp<variations; varp++) {

        cout << "init..." << endl;
        // Initialisations
        double psisq, psisq_initial;
        int randcord;
        double randno, randno2, randno3;
        double ratio;
        
        double xstep = 0.5;

        // Position matrix
        double ** x = new double*[N];
        double *xtmp= new double[6];
        for (int i=0; i<N; i++) {
            x[i] = new double[6];
        }

        // Initial position is random
        for (int j=0; j<6; j++) {
            x[0][j] = xstep * (gsl_rng_uniform(r) - 0.5);
        }

        // Calculate wave func for initial positoin
        psisq_initial = wavefuncSq(x[0], alphas[a], betas[b]);

        // Counter
        int n = 1;

        // Loop until n = N-1
        cout << "VMC..." << endl;
        //#pragma omp for private(randno, randno2, randno3, xtmp)
        //for (int n=1; n<N; n++) {
        while (n<N) {
            // Pick random coordinate [0,5] and calc ratio
            randcord = (int)floor(6*gsl_rng_uniform(r));

            // Calculate RVs
            randno = gsl_rng_uniform(r);    // in range [0,1).
            randno3= gsl_rng_uniform(r);    // in range [0,1).

            // Transform no -> [-1,1):
            randno2= 2*randno - 1.0;

            // Calculate new step:
            for (int j=0; j<6; j++) { 
                xtmp[j] = x[n-1][j];
            }
            
            // Change random coordinate and calc wave func
            xtmp[randcord] += randno2*xstep;
            psisq = wavefuncSq(xtmp, alphas[a], betas[b]);

            // Calc wavefunc ratio
            ratio = psisq / psisq_initial;

            if (ratio > randno3) {
                // Continue, 
                // Accept location:
                for (int j=0; j<6; j++) { 
                    x[n][j] = xtmp[j];
                }
                // Increment counter
                n++;
            } 
           // else {
           //     // New loc not accepted
           //     for (int j=0; j<6; j++) {
           //         x[n][j] = x[n-1][j];
           //     }
           // }
        }

        // Calculate local energy for sample
        cout << "Local energy...";
        LocalEnergyNumerical(x \
                            ,&wavefunc \
                            ,&potential \
                            ,xstep \
                            ,N \
                            ,alphas[a] \
                            ,betas[b] \
                            ,locergs[varp]);
        cout << " \t Done..." << endl;

        energies[varp] = locergs[varp][0];
        // Calculate std.deviation
        stddevs[varp] = sqrt(1./N * (locergs[varp][1] - locergs[varp][0]*locergs[varp][0]));


        cout << "Alpha " << alphas[a] << endl;
        cout << "Beta  " << betas[b]  << endl;
        cout << "Energy " << locergs[varp][0] << endl;
        cout << "energ sq " << locergs[varp][1] << endl;

        a += aplus;
        b += bplus;
        
        // Housekeeping
        delete [] x;
        delete [] xtmp;

    }

    // Output:
    //  - x axis: alphas
    //  - y axis: local energy
    //  - y axis: stddev
    double* yvals = new double[2*variations];
    #pragma omp parallel
    {
        #pragma omp for nowait
        for (int a=0; a<variations; a++) {
            yvals[a] = energies[a];
        }
        #pragma omp for nowait
        for (int b=variations; b<2*variations; b++) {
            yvals[b] = stddevs[b-variations];
        }
    }

    if (argc > 3) { 
        outputFile2(variations, 2, alphas, yvals, &argv[3]);        
    }

    // Time spent?
    t = omp_get_wtime() - t;
    cout << "Elapsed time: " << t << "s" << endl;

    // Housekeeping
    delete [] alphas;
    delete [] betas;
    delete [] energies; 
    delete [] stddevs;
    delete [] yvals;
    delete [] locergs;


    return 0;    
}


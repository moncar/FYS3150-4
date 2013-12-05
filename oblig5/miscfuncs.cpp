#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lib.h"
#include "headers.h"

using namespace std;


void RandNoGenGauss(double **randnos, double sigma, int dims, int N)
{
    // Calculate random numbers
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
    for (int i=0; i < dims; i++) {
        randnos[i] = new double[N];
    }
    cout << "Done initialising." << endl;

    // Calculate RVs using numbers from "generator" w/ "distribution".
    for (int j=0; j < dims; j++) {
        #pragma omp parallel for
        for (int i=0; i < N; i++) {
            randnos[j][i]    = gsl_ran_gaussian(r, sigma);
    //        x4[j][i]   = distribution(generator);
        }
        cout << "Done w/ round " << j << endl;
    }

}

double LocalEnergy1(double x1, double y1, double z1 \
                , double x2, double y2, double z2 \
                , double alpha, int N)
{
    // Calculates the local energy for the test wave function
    // PSI_T1
    
    double alfalfa = alpha*alpha;
    double r1sq;
    double r2sq;
    double r1r2; 
    
    for (int i=0; i < N; i++) {
        // loop over all steps

        r1sq = x1*x1 + y1*y1 + z1*z1;
        r2sq = x2*x2 + y2*y2 + z2*z2;
        r1r2 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

        //cout << r1sq+r2sq << " \t" << 1./r1r2 << endl;
        // Calculate local energy, store value in outputarray
        return 0.5 * alfalfa * (2. - alfalfa * (r1sq + r2sq)) \
              +0.5 * (r1sq + r2sq) + 1./r1r2;
    } 

}
    


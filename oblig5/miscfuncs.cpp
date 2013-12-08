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
    //omp_set_dynamic(0);
    //omp_set_num_threads(dims);
    #pragma omp parallel
    for (int j=0; j < dims; j++) {
        #pragma omp for
        for (int i=0; i < N; i++) {
            randnos[j][i]    = gsl_ran_gaussian(r, sigma);
    //        x4[j][i]   = distribution(generator);
        }
//        cout << "Done w/ round " << j << endl;
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
        return 1./r1r2;
//        return alfalfa * ( alfalfa * (r1sq + r2sq) - 6) \
//                + 0.5 * (r1sq - r2sq)   \
//                + 1./r1r2;
//        
        
//        return 0.5 * alfalfa * (2. - alfalfa * (r1sq + r2sq)) \
//              +0.5 * (r1sq + r2sq) + 1./r1r2;
    } 

}
    
double potential(double* x)
{
    // Returns potential term for wave functions
    // Input:
    //  - Array of x-positions
    // Output:
    //  - Double: potential

    double r1sq;
    double r2sq;
    double r1r2; 
    
    r1sq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    r2sq = x[3]*x[3] + x[4]*x[4] + x[5]*x[5];
    r1r2 = sqrt(  (x[0]-x[3])*(x[0]-x[3]) \
                + (x[1]-x[4])*(x[1]-x[4]) \
                + (x[2]-x[5])*(x[2]-x[5]));

    return 0.5 * (r1sq + r2sq) + 1./r1r2;
}


void LocalEnergyNumerical(double** x \
                        ,   double (*wvfunc)(double*, double, double) \
                        ,   double (*potential)(double*) \
                        ,   double xstep \
                        ,   int N   \
                        ,   double alpha \
                        ,   double beta \
                        ,   double* energies )
{
    // Calculates the local energy using the definition of the Hamiltonian
    // operator:
    //      H = -0.5 DEL^2   +   V
    //
    // The local energy is given as
    //      Eloc    =   1./psi H psi     =   -0.5/psi * DEL^2 psi + V

    double kinenergy = 0.0;
    double potenergy = 0.0;
    double totenergy = 0.0;
    double totenergysq=0.0;

    double wfunc_forw= 0.0;
    double wfunc_back= 0.0;
    double wfunc_pres= 0.0;

    // To hold copy of positions, these will be garbled
    double** xforw = new double*[N];
    double** xback = new double*[N];

    for (int i=0; i<N; i++) {
        xforw[i] = new double[6];
        xback[i] = new double[6];
    }
    for (int i=0; i<N; i++) {
        for (int j=0; j<6; j++) {
            xforw[i][j] = x[i][j];
            xback[i][j] = x[i][j];
        }
    }

    

    // Kinetic term
    for (int i=0; i<N; i++) {
        for (int j=0; j<6; j++) {
            xforw[i][j] = x[i][j] + xstep;
            xback[i][j] = x[i][j] - xstep;

            // Evaluate wavefunc:
            wfunc_back = (*wvfunc)(xback[i], alpha, beta);
            wfunc_forw = (*wvfunc)(xforw[i], alpha, beta);
            wfunc_pres = (*wvfunc)(x[i], alpha, beta);

            // Differentiate:
            kinenergy = 0.5*(wfunc_forw -2*wfunc_pres +wfunc_back);

            // Local energy: 1/psi * H psi
            kinenergy /= wfunc_pres*xstep*xstep;

            // Reset, to be reused
            xforw[i][j] = x[i][j];
            xback[i][j] = x[i][j];
        }
        potenergy = potential(x[i]);
        totenergy += (-kinenergy + potenergy);
        totenergysq+=pow((-kinenergy+potenergy), 2);
    }

    // Calculate energy average and squared
    totenergy /= N;
    totenergysq /= N;

    // Housekeeping
    delete [] xforw;
    delete [] xback;

    energies[0] = totenergy;
    energies[1] = totenergysq;

}









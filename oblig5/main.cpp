#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lib.h"
#include "headers.h"


// MHJ: Hermite zeros and weights
#define EPS 3.0e-14
#define PIM4 0.7511255444649425
#define MAXIT 10


using namespace std;





int main(int argc, char* argv[]) {

    int N;

    if (argc < 2) {
        N = 2*2;
        cout << "Launching integrator. N=" << (int)N << endl;
    }
    else {
        N = atoi(argv[1]);
        cout << "Launching integrator. N=" << (int)N << endl;
    }

    // Display number of possible threads
    int numthreads = omp_get_max_threads();
    int numcores = omp_get_num_procs();
    cout << "Number of avail threads: " << numthreads << endl;
    cout << "Number of avail cores: " << numcores << endl;
//    omp_set_dynamic(0);
//    omp_set_num_threads(numthreads);


    cout << "[1]:\tGauss-Legendre quadrature,"  << endl;
    cout << "[2]:\tGauss-Hermite quadrature,"   << endl;
    cout << "[3]:\tBrute force Monte Carlo,"    << endl;
    cout << "[4]:\tImportance sampling Monte Carlo,"<< endl;
    cout << "[5]:\tVariational Monte Carlo, wavefunc 1,"<< endl;
    int choice;
    if (argc < 3) {
        cin >> choice;
    }
    else {
        choice = atoi(argv[2]);
    }
    switch(choice) {
        case 1: 
            cout << "Entered 1" << endl;
            int statusMain1;
            statusMain1 = mainGaussLegendre(argc, argv, N);

            return statusMain1;
        case 2:
            cout << "entered 2" << endl;

            int statusMain2;
            statusMain2 = mainGaussHermite(argc, argv, N);

            return statusMain2;
        case 3:
            cout << "entered 3" << endl;

            int statusMain3;
            statusMain3 = mainMonteCarloBruteForce(argc, argv, N);

            return statusMain3;
        case 4:
            cout << "entered 4" << endl;

            int statusMain4;
            statusMain4 = mainMonteCarloImportanceSampling(argc, argv, N);

            return statusMain4;
        case 5:
            cout << "entered 5" << endl;

            int statusMain5;
            statusMain5 = mainMonteCarloVMC1(argc, argv, N);

            return statusMain4;
        default:
            cout << "Enter valid number..." << endl;
            return 1;
    }


}

#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lib.h"
#include "headers.h"

void printContents(double* xin)
{
    for (int j=0; j<6; j++) {
        cout << xin[j] << " \t";
    }

    cout << endl;
}


int mainTESTSUITE(int argc, char* argv[], int N) {
    
    cout << "Launching TEST SUITE. N=" << (int)N << endl;

    //Test LocalEnergyNumerical

    // def x vals to be determined using metropolis algo
    double** xvals = new double*[N];
    for (int i=0; i<N; i++) {
        xvals[i] = new double[6];
        for (int j=0; j<6; j++) {
            xvals[i][j] = j+1.0 + i;
        }
    }

    double xstep = 0.5;
    double alpha = 1.0;

    double* locenergy;

    printContents(xvals[0]);
    printContents(xvals[1]);
    
    LocalEnergyNumerical(xvals \
                        ,&wavefuncT1 \
                        ,&potential \
                        ,xstep \
                        , N \
                        , alpha \
                        , locenergy);

    cout << "Energy " << locenergy[0] << endl;
    cout << "energ sq " << locenergy[1] << endl;


    // Housekeeping
    delete [] xvals;

    return 0;
}

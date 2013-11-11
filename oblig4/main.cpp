#include <iostream>

#include "headers.h"

using namespace std;

template <typename T> int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}

double initcond(double x) {
    return x-1;
}


int main(int argc, char* argv[]) {

    if (argc < 4) {
        cout << "Error: Program requires four argv[i].\n \
                 argv[1]: Forward-Euler results .txt,\n \
                 argv[2]: Backward-Euler results .txt,\n \
                 argv[3]: Impl. Crank-Nicolson results .txt\n \
                 argv[4]: time-steps file .txt\n \
        \nBye-bye, " << endl;
        exit(1);
    }

    double d = 1.0;
    int    N = 10;
    double step = (d-0)/(N);
    double alpha = 0.5;
    int tsteps = 1./(alpha * step * step);
    double dt  = 1./tsteps;

    cout << tsteps << "!!!!" << endl;

    forward_euler(N, tsteps, step, alpha, &argv[1]);
    backward_euler(N, tsteps, step, alpha, &argv[2]);
    crank_nicolson(N, tsteps, step, alpha, &argv[3]);



    // Print steps to file: required to measure error
    double* XVAL = new double[N+1];
    double* TVAL = new double[tsteps];

    for (int t=0; t < tsteps; t++) {
        TVAL[t] = t * dt;
    }
    for (int i=0; i < N+1; i++) {
        XVAL[i] = i*step;
    }

    outputFile2(N+1, 1, XVAL, XVAL, &argv[4]);

    delete [] XVAL;
    delete [] TVAL;

    return 0;
}

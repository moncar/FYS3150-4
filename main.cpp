#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>

#include "lib.h"
#include "headers.h"

#define EPS 3.0e-14
#define PIM4 0.7511255444649425
#define MAXIT 10


using namespace std;


double wavefunc(double x1, double y1, double z1, \
                double x2, double y2, double z2) {
    
    double alpha = 1.0;

    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
    double r1r2=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));


    if (r1r2 > 1e-6) {
        return exp( - (alpha*alpha) * (r1*r1 + r2*r2) ) / (r1r2) ;
    }
    else {
        return 0.0;
    }
}

double wavefunc_red(double x1, double y1, double z1, \
                    double x2, double y2, double z2) {
    
    double alpha = 1.0;

    double r1r2=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

    //if (r1r2 > 1e-6) {
        return 1.0 / (r1r2+1.0) ;
    //}
    //else {
    //    return 0.0;
    //}
}

void GaussHermite(double *x,double *w, int n)
{
    int i,its,j,m;
    double p1,p2,p3,pp,z,z1;

    m=(n+1)/2;
    for (i=1;i<=m;i++) {
        if (i == 1) {
                z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
        } else if (i == 2) {
                z -= 1.14*pow((double)n,0.426)/z;
        } else if (i == 3) {
                z=1.86*z-0.86*x[1];
        } else if (i == 4) {
                z=1.91*z-0.91*x[2];
        } else {
                z=2.0*z-x[i-2];
        }
        for (its=1;its<=MAXIT;its++) {
            p1=PIM4;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
            }
            pp=sqrt((double)2*n)*p2;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "MAXIT reached!" << endl;
        x[i]=z;
        x[n+1-i] = -z;
        w[i]=2.0/(pp*pp);
        w[n+1-i]=w[i];
    }
}

int main1(int argc, char* argv[]) {
    
    int N;

    N = atoi(argv[1]);
    cout << "Launching Gauss-Legendre integrator. N=" << (int)N << endl;

    // Gauss-Legendre:  integrate wave-function
    //                  which is approx zero at r=5

    double lima = -5.0;
    double limb = +5.0;


    double* x = new double[N];  // zeros/ mesh points
    double* w = new double[N];  // weights
    double intsum = 0.0;        // integration sum

    // init: mesh points and weights
    gauleg(lima, limb, x,w, N);

    // integrate
    #pragma omp parallel for reduction(+:intsum)
    for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
    for (int k=0; k < N; k++) {
    for (int l=0; l < N; l++) {
    for (int m=0; m < N; m++) {
    for (int n=0; n < N; n++) {
        intsum += w[i]*w[j]*w[k]*w[l]*w[m]*w[n] \
                * wavefunc(x[i], x[j], x[k], x[l], x[m], x[n]);
    }}}}}}

    printf("Integral sum: %g , with N=%d \n", intsum, N);
    
    // Housekeeping
    delete [] x;
    delete [] w;
    
    
    return 0;
}

int main2(int argc, char* argv[]) {
    
    int N;

    N = atoi(argv[1]);
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


int main3(int argc, char* argv[]) {
    
    // Monte-Carlo
    
    int N;

    N = atoi(argv[1]);
    cout << "Launching Monte Carlo integrator. N=" << (int)N << endl;


    // Brute-force Monte Carlo

    // Find random variables in range: [0,1] to [-infty, +infty]
    //
    // 6 dimensions:
   
    // init random seed 
    srand(time(0));
    mt19937 eng(time(NULL));
    random_device rd;
    uniform_real_distribution<double> uniform_real(-5.0,5.0);
   
    // Jacobi determinant
    // [b-a] = 2a
    double jacobiDet  = pow((2*5), 6.0);

    double** xxx = new double*[6];
    double*  fx  = new double[N];
    #pragma omp parallel for
    for (int i=0; i < 6; i++) {
        xxx[i] = new double[N];
    }

    #pragma omp parallel for
    for (int j=0; j < 6; j++) {
        for (int i=0; i < N; i++) {
            xxx[j][i]   = uniform_real(rd);
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

    double integral = intsum3/N * jacobiDet;
    sigmasq         = sigmasq/N;
    variance        = sigmasq - integral*integral;
    variance        = variance / N ; // ???

    printf("Integral sum: %.12g , with N=%d \n", integral, N);
    cout << "\t" << sigmasq << " x^2: " << integral*integral << endl;
    printf("\tVariance: %.12g,\n\tstd.dev: %.12g\n", variance, sqrt(variance));

    if (argc > 3) { 
        outputFile2(N, 1, xxx[1], fx, &argv[3]);
    }
    
    // Housekeeping
    delete [] xxx;
    delete [] fx;

    return 0;
}






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

    cout << "[1]:\tGauss-Legendre quadrature,"  << endl;
    cout << "[2]:\tGauss-Hermite quadrature,"   << endl;
    cout << "[3]:\tBrute force Monte Carlo,"    << endl;
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
            statusMain1 = main1(argc, argv);

            return statusMain1;
        case 2:
            cout << "entered 2" << endl;

            int statusMain2;
            statusMain2 = main2(argc, argv);

            return statusMain2;
        case 3:
            cout << "entered 3" << endl;

            int statusMain3;
            statusMain3 = main3(argc, argv);

            return statusMain3;
        default:
            cout << "Enter valid number..." << endl;
            return 1;
    }


}

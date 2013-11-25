#include <iostream>
#include <math.h>
#include <omp.h>

#include "lib.h"

using namespace std;


double wavefunc(double x1, double y1, double z1, \
                double x2, double y2, double z2) {
    
    double alpha = 1.0;

    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
    double r1r2=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));


    if (r1r2 > 1e-6) {
        return exp( - (alpha*alpha) * (r1*r1 + r2*r2) ) / r1r2 ;
    }
    else {
        return 0.0;
    }
}

double wavefunc_red(double x1, double y1, double z1, \
                    double x2, double y2, double z2) {
    
    double alpha = 1.0;

    double r1r2=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));


    if (r1r2 > 1e-6) {
        return 1.0 / r1r2 ;
    }
    else {
        return 0.0;
    }
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


int main(int argc, char* argv[]) {


    int N = 2*2*2*2*2*2*2;
    double lima = -5.0;
    double limb = +5.0;

    // Gauss-Legendre:  integrate wave-function
    //                  which is approx zero at r=5

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


    // Gauss-Hermite:   integrate wave-function
    //                  using ortho func: Hermite, +/- infty
    
    double* xx = new double[N+1]; // zeros/mesh points
    double* ww = new double[N+1]; // weights
    double intsum2= 0.0;        // integration sum: G-H

    // init: mesh points and weights
    GaussHermite(x, w, N);

    // integrate

    #pragma omp parallel for reduction(+:intsum)
    for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
    for (int k=0; k < N; k++) {
    for (int l=0; l < N; l++) {
    for (int m=0; m < N; m++) {
    for (int n=0; n < N; n++) {
        intsum2 += w[i]*w[j]*w[k]*w[l]*w[m]*w[n] \
                * wavefunc_red(x[i], x[j], x[k], x[l], x[m], x[n]);
    }}}}}}

    printf("Integral sum: %g , with N=%d \n", intsum2, N);

    //



////    int treid; 
////    double sumtot = 0;
////    
////    #pragma omp parallel private(treid)
////    {
////        treid = omp_get_thread_num();
////        std::cout << "Hello World! " << treid << "\n";
////        //sumtot += wavefunc(treid, treid+1, treid, treid+2, treid, treid+3);
////        sumtot += treid;
////
////
////    }
////
////
    
    return 0;

}

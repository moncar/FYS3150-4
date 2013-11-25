#include <iostream>
#include <math.h>
#include <omp.h>

#include "lib.h"

using namespace std;


double wavefunc(double x1, double y1, double z1, \
                double x2, double y2, double z2) {
//double wavefunc(double aX[]) {
//    
//    double x1 = aX[0];
//    double y1 = aX[1];
//    double z1 = aX[2];
//    double x2 = aX[3];
//    double y2 = aX[4];
//    double z2 = aX[5];


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

int main(int argc, char* argv[]) {


    int N = 2*2*2*2*2*2;
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
    for(int n=0; n < N; n++) {
        intsum += w[i]*w[j]*w[k]*w[l]*w[m]*w[n] \
                * wavefunc(x[i], x[j], x[k], x[l], x[m], x[n]);
    }}}}}}

    printf("Integral sum: %g , with N=%d \n", intsum, N);



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

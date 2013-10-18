// Sixth order Dormand-Prince SOLVER ROUTINE

// Swallows functions, initial values and solution arrays

//#include <iostream>
//#include <cmath>
#include <armadillo>

#include "headers.h"

using namespace std;
using namespace arma;

void solver_DPRI(   arma::vec (*vFuncIn)(arma::vec, arma::vec) \
                , arma::mat *X \
                , int N \
                , int M \
                , double dStep) 
{
    // Function: vFuncIn( X(i), X(i+1) )

    // Butcher tableu
    double a11, a21, a22, a31, a32, a33, a41, a42, a43, a44;
    double a51, a52, a53, a54, a55;
    double a61,      a63, a64, a65, a66;

    // 5th order coefficients
    double b0, b2, b3, b4, b5;

    // 6th order coefficients
    double bb0, bb2, bb3, bb4, bb5, bb6;

    // k-vectors
    arma::vec k1, k2, k3, k4, k5, k6;

    a11=1./5;
    a21=3./40      ; a22 = 9./40;
    a31=44/45      ; a32 = -56./15      ; a33 = 32./9;
    a41=19372./6561; a42 = -25360./2187 ; a43 = 64448./6561;
    a44=-212./729  ;
    a51=9017./3168 ; a52 = -355./33     ; a53 = 46732./5247;
    a54=49./176    ; a55 = -5103./18656
    a61=35./384    ;                      a63 = 500./1113;
    a64=125./192   ; a65 = -2187./6784  ; a66 = 11./84;
    ## 5th order coeffs
    b0=35./384     ; 0                  ; b2 = 500./1113; b3= 125./192;
    b4= -2187./6784; b5 = 11./84;
    ## 6th rder coeffs
    bb0=5179./57600; 0                  ; bb2 =7571./16695;
    bb3=363./650   ; bb4 =-92097./339200; bb5 = 187./2100; bb6 =1./140;
 
    while i > 0:
       c1h =    float(h)/5
       c2h = 3* float(h)/10
       c3h = 4* float(h)/5
       c4h = 8* float(h)/9
       c5h =    float(h)
       c6h =    float(h)
    
       k0 = dStep* vFuncIn(ti,xi)
       k1 = dStep* vFuncIn(ti+c1h, xi + a11     *k0)
       k2 = dStep* vFuncIn(ti+c2h,xi +  a21     *k0 +a22   * k1)
       k3 = dStep* vFuncIn(ti+c3h,xi  + a31     *k0 +a32   * k1 \
                            + a33     *k2)
       k4 = dStep* vFuncIn(ti+c4h,  xi +a41     *k0 +a42   * k1 \
                            +a43      *k2 +a44   * k3)
       k5 = dStep* vFuncIn(ti+c5h, xi +a51      *k0 +a52   * k1 \
                            +a53      *k2 +a54   * k3 \
                            +a55      *k4)
       k6 = dStep* vFuncIn(ti+c6h, xi + a61     *k0 +a63   * k2 \
                            + a64     *k3 +a65   * k4     + a66 * k5)
       // print k0
       // 5th order RK:
       ya = xi + b0       *k0 + b2 * k2 \
               + b3       *k3 + b5 * k5
       // 6th order RK
       yb = xi + bb0      *k0 + bb2  * k2 \
               + bb3      *k3 + bb4  * k4 + bb5 * k5 + bb6 * k6
       yb = array(yb)
 
       J = max(shape(yb))         // Number of array elements
       errors = zeros(J,float)    // to find largest deviation
       scales = zeros(J,float)    // scale eps according to source of error
       deltas = zeros(J,float)    // calculate relative error
       for j in xrange(0,max(J, 1)): // max(J-1,1) in case ya,yb,xi: scalars
          errors[j] = abs(ya[j] - yb[j])
 	 scales[j] = eps_0*max(abs(xi[j]),abs(ya[j])) 
          deltas[j] = (float(errors[j])/scales[j])**2
       
       err1 = sqrt(1./J * sum(deltas))  // calculated error
       err0= 1.                         // desired error
       // if err <= 1: accept step, if err > 1 recalculate
       print err1
       if err1 <= 1 + 1e-8:
          // acceptable, try to increase step?
          print "INCR"
          tt = ti + h
          h = sign*min(abs(h* abs(err0 / err1)**(1./5)), abs(10*h))
          i = 0 // exit loop
       elif err1 > 1 + 1e-8:
          // unacceptable, decrease step
          print "DECR"
          h1 = h* abs(err0 / err1)**(1./5)
          print h/h1
          if float(h)/h1 >= 5:
             print float(h1)/h
             // Break loop if h grows too small
             print "EXIT"
             tt = ti + h
             h = h1
             i = 0
          else:
             h = h1
 
    //Return integrated next step: yb, next timestep: tt, appr. step size: h
    return yb,tt,h



    // Routine
    
    for (int iii=0; iii<N-1; iii++) {
        // Iterate until the end of the array
        for (int j=0; j<M; j++) {
            // Iterate through each equation
            // Recalculate each of them as dStep could change (later)
            
            arma::vec k1 = (*vFuncIn)(X->col(iii), X->col(iii+1));
            arma::vec k2 = (*vFuncIn)(X->col(iii) + k1*0.5*dStep,\
                   X->col(iii) + k1*0.5*dStep) ;
            arma::vec k3 = (*vFuncIn)(X->col(iii) + k2*0.5*dStep,\
                   X->col(iii) + k2*0.5*dStep) ;
            arma::vec k4 = (*vFuncIn)(X->col(iii) + k3*dStep, \
                   X->col(iii) + k3*dStep) ;



            X->col(iii+1) = X->col(iii) + \
                          1./6 * dStep * (k1 + 2*k2 + 2*k3 + k4);
        }
    }


}


#include <cmath>
#include "headers.h"


// Wave function declarations

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

    if (r1r2 > 1e-6) {
        return 1.0 / (r1r2) ;
    }
    else {
        return 0.0;
    }
}


double wavefuncT1(double* x, double alpha, double) 
{
    
    double r1sq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r2sq = x[3]*x[3] + x[4]*x[4] + x[5]*x[5];

    return exp(- alpha*alpha * (r1sq + r2sq )/2.0);
}

double wavefuncT1sq(double* x, double alpha, double) 
{
    
    double r1sq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r2sq = x[3]*x[3] + x[4]*x[4] + x[5]*x[5];

    return exp(- alpha*alpha * (r1sq + r2sq ));
}

double wavefuncT2(double* x, double alpha, double beta)
{
    double r1sq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r2sq = x[3]*x[3] + x[4]*x[4] + x[5]*x[5];
    double r1r2 = sqrt((x[0]-x[3])*(x[0]-x[3]) + (x[1]-x[4])*(x[1]-x[4]) + (x[2]-x[5])*(x[2]-x[5]));

    return exp(- alpha*alpha * (r1sq + r2sq )/2.0 \
               + r1r2/(2.0 * (1.0 + beta* r1r2) ) );
}

double wavefuncT2sq(double* x, double alpha, double beta)
{
    double r1sq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r2sq = x[3]*x[3] + x[4]*x[4] + x[5]*x[5];
    double r1r2 = sqrt((x[0]-x[3])*(x[0]-x[3]) + (x[1]-x[4])*(x[1]-x[4]) + (x[2]-x[5])*(x[2]-x[5]));

    return exp(- alpha*alpha * (r1sq + r2sq ) \
               + r1r2/(1.0 + beta* r1r2) );
}


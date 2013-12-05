#include <iostream>
#include <math.h>

#include "lib.h"
#include "headers.h"


// MHJ: Hermite zeros and weights
#define EPS 3.0e-14
#define PIM4 0.7511255444649425
#define MAXIT 10


using namespace std;

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

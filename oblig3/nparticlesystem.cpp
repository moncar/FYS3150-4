// N-PARTICLE SYSTEM FORCE/ACCELERATION EXPRESSION

#include <iostream>
#include <armadillo>

#include "headers.h"

using namespace std;
using namespace arma;

const double pi = 4*atan(1);

arma::vec aSystem(double dTi, arma::vec Xi) {
    // Set up vector containing velocities, positions and time
    // for entire system
    // Pos: (x,y) -> 2 * nPlanets
    // Vel: (x,y) -> 2 * nPlanets
    // Mass:      -> 1 * nPlanets
    // time          1
    
    //int nPlanets = (sizeof(Xi)/sizeof(double) - 1) / (5);
    double nPlanets = ( Xi.size() )/5;
    vec Y(5*nPlanets);

    //for planet k
    //    interacting with planet j != k:
    

    for (int k=0; k<nPlanets; k++) {
        //double fFtot = 0;   // Total force on each planet k
        double fFtotx= 0;   //  x-component
        double fFtoty= 0;   //  y-component
        double fFloc = 0;   // Force on planet k from planet j
        double fFlocx= 0;   //  x-component
        double fFlocy= 0;   //  y-component
        
        for (int j=0; j<nPlanets; j++) {
            if (j != k) {
                // Difference in position
                double fDX = Xi[j*5+2] - Xi[k*5+2];
                double fDY = Xi[j*5+3] - Xi[k*5+3];
                double fRad= sqrt(fDX*fDX + fDY*fDY);

                // Forces, divided by M[k]
                // Newtonian:
                fFloc = Xi[j*5+4] / (fRad * fRad);

                fFlocx= fFloc * (Xi[j*5+2] - Xi[k*5+2])/fRad;  
                fFlocy= fFloc * (Xi[j*5+3] - Xi[k*5+3])/fRad;
                
                fFtotx += fFlocx;
                fFtoty += fFlocy;

            }
        }
        // d2x / dt2
        // Fx/M[k] = F/M[k] * cos theta = F/M[k] * x/r 
        Y[k*5]      = fFtotx; 
        
        // d2y / dt2
        // Fy/M[k] = F/M[k] * cos theta = F/M[k] * y/r 
        Y[k*5+1]    = fFtoty;

        // dx  / dt
        Y[k*5+2]    = Xi[k*5] ;
        // dy  / dt
        Y[k*5+3]    = Xi[k*5+1] ;

        // dm / dt
        Y[k*5+4]    = 0;
    }

    return Y;
}

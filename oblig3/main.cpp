#include <iostream>
#include <cmath>
#include <armadillo>

#include "headers.h"
#include "SolarSystem.hpp"

using namespace std;
using namespace arma;

template <typename T> int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}


inline double fVel(double fMass, double fRadius) {
    // Calculates the velocity given by the expression
    // for the centripetal acceleration
    return sqrt(fMass/fRadius);
}

int main(int argc, char* argv[]) {
    
    int N= 3650*1000; // 414 Mercury years 414*10*88; //365/2;
    //int N=5;
    double pi = 4*atan(1);
    double dt = 2*pi/365/10;
    double t0 = 0.0;

    SolarSystem Sol;
    double fMsol = 332948.6; // Earth masses


    // // //    SETTING UP SYSTEM   // // //
    
    // Adding celestial objects using:
    // Sol.AddPlanet(Mass,      xpos,   ypos,   xvel,   yvel)
   
    // Sun
    //Sol.AddPlanet(1.,           0.,       0.,   0.,     0.); 
    // Mercury
    //Sol.AddPlanet(0.055/fMsol,  0.387098, 0.,   0., -fVel(1., 0.387098));
    // Mercury at perihelion, with given velocity:
    Sol.AddPlanet(0.055/fMsol,  0.3075,   0.,   0.,  12.44/(2*pi));
    // Venus
    Sol.AddPlanet(0.815/fMsol, -0.723327, 0.,   0., -fVel(1., 0.723327));
    // Earth
    Sol.AddPlanet(1./fMsol,     1.,     0.,     0.,  fVel(1., 1.));
    // Moon
    Sol.AddPlanet(0.0123/fMsol, 1.00257,0.,0., 1.+ fVel(1./fMsol, 0.00257));
    // Mars
    Sol.AddPlanet(0.107/fMsol, -1.523679, 0., 0., fVel(1., 1.523679));
    // Jupiter
    Sol.AddPlanet(95.152/fMsol, 0., 9.582, fVel(1., 9.582), 0.);
    // Saturn
    Sol.AddPlanet(317.8/fMsol, 0., 5.204267, fVel(1., 5.204267), 0.);
    // Uranus
    Sol.AddPlanet(14.535/fMsol, 0., 19.229, -fVel(1., 19.229), 0.);
    // Neptune
    Sol.AddPlanet(17.147/fMsol, 0., -30.103, fVel(1., 30.103), 0.);
    
    // We are missing the Sun!
    //  -> Find total momentum
    //  -> Determine initial conditions of Sun opposite to rest of System
    
    // Find momentum of system:
    double Px, Py, Mtot;
    Sol.TotMomentum(Px, Py, Mtot);
    cout << Px << " " << Py << " " << Mtot << endl;

    // Require that Sun has tot momentum opposite to system
    double VxSun, VySun, mSun;
    mSun = 1.;
    VxSun = -1./mSun * Px;
    VySun = -1./mSun * Py;

    // Define Sun:
    Sol.AddPlanet(mSun,         0.,     0.,     VxSun,  VySun);

    // // //    DONE: SETTING UP SYSTEM     // // //


    // Set up matrix A: holds initial values for all planets
    mat A = Sol.ConstructArray(N);

    cout << "Class has # planets: " << Sol.PlanetCounter() << endl;

    // Send matrix to solver, the class sends the general n-body function
    // defined in NPARTICLESYSTEM.CPP
    Sol.SolveAll(A, dt, t0);

    int nPlanets = Sol.PlanetCounter();    

    //A.print();

    
    // // //    SAVE TO FILE    // // //
    // Pass time, position to file
    
    double* POS = new double[2*nPlanets*N];
    double* TIME= new double[N];

    for (int iii=0; iii<N; iii++) {
        TIME[iii] = t0 + iii*dt;
    }

    for (int ppp=0; ppp<nPlanets; ppp++) {
        for (int kkk=0; kkk<N; kkk++) {
            int subscript = 2*ppp;
            POS[kkk +  subscript *N] = A(ppp*5+2,kkk);
            POS[kkk +  (subscript+1) *N] = A(ppp*5+3,kkk);
            
            //TIME[kkk]   = A(5*nPlanets,kkk);
        }
    }

    if (argc < 2) {
        cout << "ERROR! Missing output file specification:\n \
                \t ./obl3.x outputfile.txt" << endl;
        exit(1);
    }
    else {
        outputFile(N, 2*nPlanets, TIME, POS, &argv[1]);
    }
    
    // // //    DONE: SAVE TO FILE    // // //

    // Housekeeping
    delete [] POS;
    delete [] TIME;


    return 0;
}

#include <iostream>
#include <armadillo>

#include "headers.h"
#include "SolarSystem.hpp"

using namespace std;
using namespace arma;

void SolarSystem::AddPlanet(double mass=1, double xpos=1, double ypos=0 \
    ,double xvel=0, double yvel=0)
{
    nPlanets += 1;
    M.push_back(mass);
    vx.push_back(xvel);
    vy.push_back(yvel);
    x.push_back(xpos);
    y.push_back(ypos);
    
}

arma::mat SolarSystem::ConstructArray(int N)
// Set up matrix to hold all velocities and positions
{
    nIterations = N;
    arma::mat A(5*nPlanets, N);
    
    // Data:
    // Columns: 0,1: xvel, yvel -> 2*nPlanets
    //          2,3: xpos, ypos -> 2*nPlanets
    //            4: time       -> 1
    //
    // Each planet is assigned an index k
    
    // Fill in: initial values

    for (int k=0; k<nPlanets; k++) {
        A(k*5  ,0) = vx[k];
        A(k*5+1,0) = vy[k];
        A(k*5+2,0) =  x[k];
        A(k*5+3,0) =  y[k];
        A(k*5+4,0) =  M[k];
    }
    A(nPlanets, 0) =    0;  // init. time

    return A;
}

inline double SolarSystem::Radius(double x, double y) {
    return sqrt(x*x + y*y);
}

void SolarSystem::TotMomentum(double& dPX, double& dPY, double& dM)
    // Return total momentum for composite system from initial conditions
    // Arguments will be changed
    // dPX = COM momentum in x-dir
    // dPY = COM momentum in y-dir
    // dM  = total mass
{
    // Find momentum:
    // P = sum_i (m_i v_i)
    for (int i=0; i<nPlanets; i++) {
        dPX += M[i] * vx[i];
        dPY += M[i] * vy[i];
        dM  += M[i];
    }
}

void SolarSystem::PotEnergy(arma::mat &Y     /* Holding all pos&vel  */ \
                        ,   arma::vec &aCOMx /* To hold: COM x-coord */ \
                        ,   arma::vec &aCOMy /* To hold: COM y-coord */ \
                        ,   arma::vec &aPot  /* To hold: calc pot. energy*/\
                        ,   int nIterations  /* No of iterations */\
                        ,   int nPlanetNo    /* Chosen planet/object */ )
    // Find potential energy at all time steps for planet nPlanetNo
    // through a calculation of:
    //      - COM for all other planets
    //      - Radius: COM <-> planet nPlanetNo
    //      - Finding potential energy between two bodies
{
    
    // Find total mass for system, except planet to be analysed:
    for (int i=0; i<nPlanets; i++) {
       fMtot += M[i];
    }
    fMtot -= M[nPlanetNo];

    // Calculate for each time step:
    for (int i=0; i<nIterations; i++) {
        COMx = 0.;
        COMy = 0.;
        for (int j=0; j<nPlanets; j++) {
            COMx += M[j] * Y(j*5+2, i);
            COMy += M[j] * Y(j*5+3, i); 
        }
        // Subtract contribution from the planet to be analysed:
        COMx -= M[nPlanetNo] * Y(nPlanetNo*5+2, i);
        COMy -= M[nPlanetNo] * Y(nPlanetNo*5+3, i);
        
        // Assign COM 
        aCOMx(i) = 1./fMtot * COMx;
        aCOMy(i) = 1./fMtot * COMy;

        // Radius: COM <-> planet nPlanetNo
        tmpRadCOM_planet    = sqrt(aCOMx(i)*aCOMx(i) + aCOMy(i)*aCOMy(i) ) \
                            - sqrt(Y(nPlanetNo*5+2,i) * Y(nPlanetNo*5+2,i)\
                                +  Y(nPlanetNo*5+3,i) * Y(nPlanetNo*5+3,i));

        // Finding potential energy between two bodies
        aPot(i) = fMtot * M[nPlanetNo] / tmpRadCOM_planet; 
    }
    
}



void SolarSystem::SolveAll(arma::mat &Y, double dt, double dt0) {
    //solver_EC(&aSystem, &Y, nIterations, 5*nPlanets, dt);
    solver_RK4(&aSystem, &Y, nIterations, 5*nPlanets, dt0, dt);
    //solver_DPRI(&aSystem, &Y, nIterations, 5*nPlanets, dt0, dt);
}


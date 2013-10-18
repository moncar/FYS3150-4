#include <iostream>
#include <cmath>
#include <armadillo>

#include "headers.h"

using namespace std;
using namespace arma;

//inline arma::vec vF1(arma::vec X) {
//   return x+2.0;
//}

template <typename T> int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}

class SolarSystem
{
private:
    // Initial values have to be stored in a list
    vector<double> M ;
    vector<double> vx;
    vector<double> vy;
    vector<double> x ;
    vector<double> y ;

    int nPlanets = 0;   // will change as heavenly objects are added
    int nIterations=0;  // to be defined by ConstructArray

public:
    void AddPlanet(double mass=1, double xpos=1, double ypos=0 \
        ,double xvel=0, double yvel=0)
    {
        nPlanets += 1;
        M.push_back(mass);
        vx.push_back(xvel);
        vy.push_back(yvel);
        x.push_back(xpos);
        y.push_back(ypos);
        
    }

    int PlanetCounter() {
        int nPlanetsPublic = nPlanets;
        return nPlanetsPublic;
    }

    arma::mat ConstructArray(int N)
    // Set up matrix to hold all velocities and positions
    {
        nIterations = N;
        arma::mat A(5*nPlanets+1, N);
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

    inline double Radius(double x, double y) {
        return sqrt(x*x + y*y);
    }
    

    void SolveAll(arma::mat &Y, double dt) {
        //solver_EC(&aSystem, &Y, nIterations, 5*nPlanets+1, dt);
        solver_RK4(&aSystem, &Y, nIterations, 5*nPlanets+1, dt);
    }


};

inline double fVel(double fMass, double fRadius) {
    return sqrt(fMass/fRadius);
}

arma::vec aSystem(arma::vec Xi, arma::vec Xii) {
    // Set up vector containing velocities, positions and time
    // for entire system
    // Pos: (x,y) -> 2 * nPlanets
    // Vel: (x,y) -> 2 * nPlanets
    // Mass:      -> 1 * nPlanets
    // time          1
    
    //int nPlanets = (sizeof(Xi)/sizeof(double) - 1) / (5);
    double nPlanets = ( Xi.size() - 1)/5;
    vec Y(5*nPlanets+1);

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
        //Y[k*5+1]    = fFtot * (Xi[j*5+3] - Xi[k*5+3])/fRad;
        Y[k*5+1]    = fFtoty;

        // dx  / dt
        Y[k*5+2]    = Xii[k*5] ;
        // dy  / dt
        Y[k*5+3]    = Xii[k*5+1] ;

        // dm / dt
        Y[k*5+4]    = 0;
    }

    // dt / dt
    Y[nPlanets*5] = 1.;

    return Y;
}

int main(int argc, char* argv[]) {
    
    int N=100000;
    double dt = 1./3650;

    SolarSystem Sol;
    double fMsol = 332946.; // Earth masses

    // Adding celestial objects using:
    // Sol.AddPlanet(Mass, xpos, ypos, xvel, yvel)
    
    // Sun
    Sol.AddPlanet(1.,           0.,     0.,     0.,       0.);
    // Mercury
    Sol.AddPlanet(0.055/fMsol,  0.387098, 0.,   0., -fVel(1., 0.387098));
    // Venus
    Sol.AddPlanet(0.815/fMsol, -0.723327, 0.,   0., -fVel(1., 0.723327));
    // Earth
    Sol.AddPlanet(1./fMsol,     1.,     0.,     0.,  fVel(1., 1.));
    // Moon
    //Sol.AddPlanet(0.0123/fMsol, 1.00257, 0., 0., fVel(1., 1.00257));
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
    mat A = Sol.ConstructArray(N);

    cout << "Class has # planets: " << Sol.PlanetCounter() << endl;

    Sol.SolveAll(A, dt);

    int nPlanets = Sol.PlanetCounter();    

    //A.print();

    
    // SAVE TO FILE
    // Pass time, position to file
    
    double* POS = new double[2*nPlanets*N];
    double* TIME= new double[N];

    for (int ppp=0; ppp<nPlanets; ppp++) {
        for (int kkk=0; kkk<N; kkk++) {
            int subscript = 2*ppp;
            POS[kkk +  subscript *N] = A(ppp*5+2,kkk);
            POS[kkk +  (subscript+1) *N] = A(ppp*5+3,kkk);
            
            TIME[kkk]   = A(5*nPlanets,kkk);
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

    // Housekeeping
    delete [] POS;
    delete [] TIME;


    return 0;
}

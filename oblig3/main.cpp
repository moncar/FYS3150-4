#include <iostream>
#include <list>
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

    vector<double>::iterator Mit  = M.begin();
    vector<double>::iterator vxit = vx.begin();
    vector<double>::iterator vyit = vy.begin();
    vector<double>::iterator xit  = x.begin();
    vector<double>::iterator yit  = y.begin();

    int nPlanets = 0;   // will change as heavenly objects are added
    int nIterations=0;  // to be defined by ConstructArray

public:
    void AddPlanet(double mass=1, double xpos=1, double ypos=0 \
        ,double xvel=0, double yvel=0)
    {
        nPlanets += 1;
        M.insert    (Mit,  mass);
        vx.insert   (vxit, xvel);
        vy.insert   (vyit, yvel);
        x.insert    (xit,  xpos);
        y.insert    (yit,  ypos);
        
        Mit++; vxit++; vyit++; xit++; yit++;
    }

    arma::mat ConstructArray(int N)
    // Set up matrix to hold all velocities and positions
    {
        nIterations = N;
        arma::mat A(4*nPlanets+1, N);
        // Data:
        // Columns: 0,1: xvel, yvel -> 2*nPlanets
        //          2,3: xpos, ypos -> 2*nPlanets
        //            4: time       -> 1
        //
        // Each planet is assigned an index k
        
        // Fill in: initial values

        for (int k=0; k<nPlanets; k++) {
            A(k*4  ,0) = vx[k];
            A(k*4+1,0) = vy[k];
            A(k*4+2,0) =  x[k];
            A(k*4+3,0) =  y[k];
        }
        A(nPlanets, 0) =    0;  // init. time

        return A;
    }

    inline double Radius(double x, double y) {
        return sqrt(x*x + y*y);
    }
    
    arma::vec aSystem(arma::vec Xi, arma::vec Xii) {
        // Set up vector containing velocities, positions and time
        // for entire system
        // Pos: (x,y) -> 2 * nPlanets
        // Vel: (x,y) -> 2 * nPlanets
        // time          1
        vec Y(4*nPlanets+1);

        //for planet k
        //    interacting with planet j != k:

        for (int k=0; k<nPlanets; k++) {
            for (int j=0; j<nPlanets; j++) {
                if (j != k) {
                    // Difference in position
                    double fDX = Xi[j*4+2] - Xi[k*4+2];
                    double fDY = Xi[j*4+3] - Xi[k*4+3];
                    double fRad= sqrt(fDX*fDX + fDY*fDY);

                    // Force direction
                    double sgnX = sgn(fDX);
                    double sgnY = sgn(fDY);

                    // Forces, divided by M[k]
                    double fFtot = M[j] / (fRad * fRad);
                    //double fFtot = 1./ fR*fR * sgn(0. - 1.) ;
                    
                    // d2x / dt2
                    // Fx/M[k] = F/M[k] * cos theta = F/M[k] * x/r 
                    Y[k*4]      = sgnX * fFtot * Xi[k*4+2]/fRad;  
                    
                    // d2y / dt2
                    // Fy/M[k] = F/M[k] * cos theta = F/M[k] * y/r 
                    Y[k*4+1]    = sgnY * fFtot * Xi[k*4+3]/fRad;

                    // dx  / dt
                    Y[k*4+2]    = Xii[k*4] ;
                    // dy  / dt
                    Y[k*4+3]    = Xii[k*4+1] ;
                }
            }
        }

        // dt
        Y[nPlanets*4] = 1.;

        return Y;
    }

    void SolveAll(arma::mat Y, double dt) {
        solver_EC(&SolarSystem::aSystem, &Y, nIterations, 4*nPlanets+1, dt);
    }


};

inline double fVel(double fMass, double fRadius) {
    return sqrt(fMass/fRadius);
}


int main(int argc, char* argv[]) {

    SolarSystem Sol;
    Sol.AddPlanet();
    Sol.ConstructArray(10);

    int N=1200;
    int M=5;
    double dt = 0.02;

    // Matrix to hold vectors:
    // M rows:  0,1: Velocity
    //          2,3: Position
    //            4: Time
    // N columns : updated for each timestep dt

    mat Y(M,N);
    Y.zeros();

    // Scales, EARTH AND SUN
    double fMySol      = 1.; // [solar masses]
    double fChiEarth   = 1.; // [AU]

    // Initial conditions, EARTH AND SUN
    Y(0,0) = 0.;   // vel
    Y(1,0) = fVel(fMySol, fChiEarth);   // vel
    Y(2,0) = 1.;                 // pos
    Y(3,0) = 0.;                 // pos
    Y(4,0) = 0.;                        // time

    //Y.print();
    //solver_EC(&aSystem, &Y, N, M, dt);

    //Y.print();



    // Pass time, position to file
    
    double* POS = new double[2*N];
    double* TIME= new double[N];

    for (int kkk=0; kkk<N; kkk++) {
        POS[kkk]   = Y(2,kkk);
        POS[kkk+N] = Y(3,kkk);
        TIME[kkk]   = Y(4,kkk);
    }

    if (argc < 2) {
        cout << "ERROR! Missing output file specification:\n \
                \t ./obl3.x outputfile.txt" << endl;
        exit(1);
    }
    else {
        outputFile(N, 2, TIME, POS, &argv[1]);
    }

    // Housekeeping
    delete [] POS;
    delete [] TIME;


    return 0;
}

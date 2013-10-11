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

    int nPlanets = 0;

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

    arma::cube ConstructArray(int N)
    // Set up cube to hold all velocities and positions
    {
        arma::cube A(6, N, nPlanets);
        // Data:
        // Columns: 0,1: xvel, yvel
        //          2,3: xpos, ypos
        //            4: time
        //            5: mass
        
        // Fill in: initial values
        for (int k=0; k<nPlanets; k++) {
            A(0,0,k) = vx[k];
            A(1,0,k) = vy[k];
            A(2,0,k) =  x[k];
            A(3,0,k) =  y[k];
            A(4,0,k) = 0.;
            A(5,0,k) =  M[k];
        }

        return A;
    }

    inline double Radius(double x, double y) {
        return sqrt(x*x + y*y);
    }

    void Paths(arma::cube* A)
    // Calculate velocities and paths
    {
        // Forces:
        double* fFk = new double[nPlanets];
        for (int k=0; k<nPlanets; k++) {
            fFk[k] = 0.0; 
        }

        for (int k=0; k<nPlanets; k++) {
            // Iterate through each src planet
            for (int j=0; j<nPlanets; j++) {
                // Iterate through the other planets
                if (k != j) {
                    fFk[k] = A(0,0,j) / \
                    ( Radius(A(2,i,k),A(3,i,k)) *  \
                      Radius(A(2,i,k),A(3,i,k)) ) * \
                    sgn(Radius(A(2,i,k),A(3,i,k)) - \
                        Radius(A(2,i,k),A(3,i,k) ));
                }
            }
        } 


    }

};

inline arma::vec aSystem(arma::vec Xi, arma::vec Xii) {
    vec Y(5);

    // Radius
    double fR = sqrt( Xi[2]*Xi[2] + Xi[3]*Xi[3] ) ;

    // Forces
    double fFtot = 1./ fR*fR * sgn(0. - 1.) ;
    // d2x / dt2
    Y[0] =  fFtot / 1. * Xi[2]/fR;      // Fx = F * cos theta = F * x/r 
    // d2y / dt2
    Y[1] =  fFtot / 1. * Xi[3]/fR;      // Fy = F * sin theta = F * y/r 

    // dx  / dt
    Y[2] = Xii[0] ;
    // dy  / dt
    Y[3] = Xii[1] ;

    // dt
    Y[4] = 1.;

    return Y;
}

inline double fVel(double fM, double fX) {
    return sqrt(fM/fX);
}

int main(int argc, char* argv[]) {

    SolarSystem Sol;
    Sol.AddPlanet();
    Sol.ConstructArray(10);

    int N=1200000;
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
    solver_EC(&aSystem, &Y, N, M, dt);

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


#include <armadillo>

#include "headers.h"

using namespace arma;

//inline arma::vec vF1(arma::vec X) {
//   return x+2.0;
//}

template <typename T> int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}

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


    return 0;
}

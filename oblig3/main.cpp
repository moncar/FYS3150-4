
#include <armadillo>

#include "headers.h"

using namespace arma;

//inline arma::vec vF1(arma::vec X) {
//   return x+2.0;
//}

inline arma::vec aSystem(arma::vec Xi, arma::vec Xii) {
    vec Y(3);
    // d2x / dt2
    Y[0] = Xi[1];  // 1./ (Xi[1] * Xi[1]) ;

    // dx  / dt
    Y[1] = Xii[0] ;

    // dt
    Y[2] = 1.;

    return Y;
}

int main(int argc, char* argv[]) {

    mat Y(3,6);
    Y.zeros();

    int N=6;
    int M=3;
    double dt = 0.2;

    // Initial conditions
    Y(0,0) = 0.;    // vel
    Y(1,0) = 1.;    // pos
    Y(2,0) = 0.;    // time

    Y.print();
    solver_EC(&aSystem, &Y, N, M, dt);

    Y.print();

    return 0;
}

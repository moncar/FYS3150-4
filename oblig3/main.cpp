
#include <armadillo>

#include "headers.h"

using namespace arma;

inline arma::vec vF1(arma::vec X) {
    return X+2.0;
}


int main(int argc, char* argv[]) {

    vec Y(3);
    Y.fill(2.0);

    vec Z(3);
    Z.zeros();

    int N=3;

    solver_EC(&vF1, &Y, &Z, N);

    return 0;
}

#include <iostream>
#include <cmath>
#include <armadillo>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>

using namespace std;
using namespace arma;

template <typename T> int sgn(T val) {
            return (T(0) < val) - (val < T(0));
}

double initcond(double x) {
    return x-1;
}

void backward_euler(int N, int tsteps, double delta_x, double alpha)
{
    // // //
    // Implicit scheme
    // // //

    double a, b, c;     // holding constants: diag;  up diag;  low diag
    
    gsl_vector* diag        = gsl_vector_calloc(N+1);
    gsl_vector* abovediag   = gsl_vector_calloc(N);
    gsl_vector* belowdiag   = gsl_vector_calloc(N);
    gsl_vector* rhs         = gsl_vector_calloc(N+1);
    gsl_vector* solution    = gsl_vector_calloc(N+1);

    // Diagonal matrix: elements
    a = c = -alpha;
    b = 1 + 2*alpha;

    // Initial conditions
    for (int i=1; i < N; i++) {
        gsl_vector_set( rhs, i, initcond(delta_x * i) );
    }
    for (int i=0; i < N; i++) { 
        // fill values: above and below diag
        gsl_vector_set( abovediag, i, a);
        gsl_vector_set( belowdiag, i, c);
    }
    // Fill diag. values
    for (int i=0; i <= N; i++) {
        gsl_vector_set( diag, i, b);
    }

    // Boundary conditions
    gsl_vector_set( rhs, N, 0.0);
    gsl_vector_set( rhs, 0, 0.0);
    gsl_vector_set( solution, 0, 0.0);
    gsl_vector_set( solution, N, 0.0);

    printf("abovediagN-1 = %g\n", gsl_vector_get( belowdiag, 0));

    // Time iteration
    for (int t=1; t<=tsteps; t++) {
        // solve set of tridiag eqs
        gsl_linalg_solve_tridiag(diag, abovediag, belowdiag, rhs, solution);

        // boundary conditions
        gsl_vector_set(solution, 0, 0.0);
        gsl_vector_set(solution, N, 0.0);

        // replace old time-solution with new
        for (int i = 0; i<= N; i++) {
            gsl_vector_set( rhs, i,  gsl_vector_get( solution, i) );
        }
    }
    
    // print solution?
    for (int i=0; i<=N; i++) {
        printf("u_%d = %g\n", i,  gsl_vector_get( solution, i));
    }


    // HOUSEKEEPING
    gsl_vector_free(diag);
    gsl_vector_free(abovediag);
    gsl_vector_free(belowdiag);
    gsl_vector_free(rhs);
    gsl_vector_free(solution);

}

void crank_nicolson(int N, int tsteps, double delta_x, double alpha) {
    // // //
    // Implicit scheme: Crank-Nicolson
    // // //
    

    // Initialise
    double a,b,c;   // a, c: off-diagonal
                    //  b  : diagonal
    
    gsl_vector* diag        = gsl_vector_calloc(N+1);
    gsl_vector* abovediag   = gsl_vector_calloc(N);
    gsl_vector* belowdiag   = gsl_vector_calloc(N);
    gsl_vector* rhs         = gsl_vector_calloc(N+1);
    gsl_vector* r           = gsl_vector_calloc(N+1);
    gsl_vector* solution    = gsl_vector_calloc(N+1);

    // Diagonal matrix: elements
    a = c = -alpha;
    b = 1 + 2*alpha;

    for (int i=0; i < N; i++) { 
        // fill values: above and below diag
        gsl_vector_set( abovediag, i, a);
        gsl_vector_set( belowdiag, i, c);
    }
    for (int i=0; i <= N; i++) {
        // Fill diag. values
        gsl_vector_set( diag, i, b);
    }


    //// // // EXPLICIT PART:
    //// // // Calculate RHS
    //gsl_vector* unew    = gsl_vector_calloc(N+1); // Next timestep

    //// Boundary conditions:
    ////      u(0,t) = 1 , t>0
    ////      u(d,t) = 0 , t>0

    //gsl_vector_set( solution, 0, 0.0);
    //gsl_vector_set( solution, N, 0.0);

    ////  Initialisation
    //for (int i=1; i < N; i++) {
    //    // init condition
    //    gsl_vector_set( solution, i, initcond( i*delta_x ) );
    //    // init new vector
    //    gsl_vector_set( unew, i, 0.0);
    //}

    //// time integration
    //for (int t = 1; t <= tsteps; t++) {
    //    for (int i = 1; i < N; i++) {
    //        // diff eq
    //        gsl_vector_set( unew, i, \
    //                  alpha      * gsl_vector_get( solution, i-1) \
    //                + (1-2*alpha)* gsl_vector_get( solution, i ) \
    //                + alpha      * gsl_vector_get( solution, i+1) \
    //                );
    //    }
    //    solution = unew;
    //}

    //// // //

    // Boundary conditions
    gsl_vector_set( rhs, N, 0.0);
    gsl_vector_set( rhs, 0, 1.0);
    gsl_vector_set( solution, 0, 0.0);
    gsl_vector_set( solution, N, 0.0);

    printf("abovediagN-1 = %g\n", gsl_vector_get( belowdiag, 0));


    // Time iteration
    for (int t = 1; t<= tsteps; t++) {
        for (int i = 1; i< N; i++) {
            gsl_vector_set( rhs, i, \
                  alpha         *gsl_vector_get( solution, i-1) \
                + (2 - 2*alpha) *gsl_vector_get( solution, i) \
                + alpha         *gsl_vector_get( solution, i+1) \
                );
            }

        gsl_vector_set( rhs, 0, 0.0);
        gsl_vector_set( rhs, N, 0.0);

        // Solve tridiag matrix
        gsl_linalg_solve_tridiag(diag, abovediag, belowdiag, rhs, solution);

        gsl_vector_set( solution, 0, 0.0);
        gsl_vector_set( solution, N, 0.0);
        
        }

    // print solution?
    printf("CRANK-NICOLSON SCHEME\n");
    for (int i=0; i<=N; i++) {
        printf("u_%d = %g\n", i,  gsl_vector_get( solution, i));
    }


    // HOUSEKEEPING
    gsl_vector_free(diag);
    gsl_vector_free(abovediag);
    gsl_vector_free(belowdiag);
    gsl_vector_free(rhs);
    gsl_vector_free(r);
    gsl_vector_free(solution);
}

int main(int argc, char* argv[]) {

    // // //
    // Explicit scheme
    // // //

    double d = 1.0;
    double N = 10;
    double step = (d-0)/(N+1);
    double alpha = 0.5;
    double tsteps = 1./(alpha * step * step);


    vec u(N+1);         // Holding current  timestep
    vec unew(N+1);      // Holding next     timestep

    // Boundary conditions:
    //      u(0,t) = 1 , t>0
    //      u(d,t) = 0 , t>0

    u(0) = unew(0) = 1.0;
    u(N) = unew(N) = 0.0;

    //  Initialisation
    for (int i=1; i < N; i++) {
        double x = i*step;
        // init condition
        u(i) = initcond(x);
        // init new vector
        unew(i) = 0;
    }

    // time integration
    for (int t = 1; t <= tsteps; t++) {
        for (int i = 1; i < N; i++) {
            // diff eq
            unew(i) = alpha*u(i-1) + (1 - 2*alpha)*u(i) + alpha*u(i+1);
        }
        u = unew;
    }




    cout << "Time steps: "<< tsteps << " with length " << 1./tsteps << " and dx = " << step << endl;
    unew.print();
    
    //backward_euler(N, tsteps, step, alpha);
    crank_nicolson(N, tsteps, step, alpha);

    return 0;
}

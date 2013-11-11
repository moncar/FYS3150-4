
#include <iostream>
#include <cmath>
#include <armadillo>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "headers.h"

using namespace std;
using namespace arma;

void forward_euler(int N, int tsteps, double delta_x, double alpha, char* filename[])
{
    // // //
    // Forward Euler
    // // //


    vec u(N+1);         // Holding current  timestep
    vec unew(N+1);      // Holding next     timestep
    mat uall(N+1, tsteps); // Holding all   timesteps

    // Boundary conditions:
    //      u(0,t) = 1 , t>0
    //      u(d,t) = 0 , t>0

    u(0) = unew(0) = 0.0;
    u(N) = unew(N) = 0.0;

    //  Initialisation
    for (int i=1; i < N; i++) {
        double x = i*delta_x;
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
        uall.col(t-1) = unew;
    }

    // // Save to file
    // 
    // Need monster-arrays with: 1) x-values, 2) computed u

    double* XVALS = new double[N+1];
    double* COMPU = new double[tsteps*(N+1)];

    for (int i=0; i<=N; i++) {
        // iterate through each x-element
        XVALS[i] = i*delta_x;
    }
    for (int t=0; t<tsteps; t++) {
        // iterate through each timestep
        for (int i=0; i<=N; i++) {
            COMPU[i + t*(N+1)] = uall(i, t);
        }
    }

    outputFile2(N+1, tsteps, XVALS, COMPU, &filename[0]);

    delete [] XVALS;
    delete [] COMPU;    

    // // // Done saving


    cout << "Time steps: "<< tsteps << " with length " << 1./tsteps << " and dx = " << delta_x << endl;
    unew.print();

}

void backward_euler(int N, int tsteps, double delta_x, double alpha, char* filename[])
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
    gsl_matrix* uallIMPL      = gsl_matrix_alloc(tsteps,N+1);

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
    gsl_matrix_set_row(uallIMPL, 0, solution);

    printf("abovediagN-1 = %g\n", gsl_vector_get( belowdiag, 0));

    // Time iteration
    for (int t=1; t<tsteps; t++) {
        // solve set of tridiag eqs
        gsl_linalg_solve_tridiag(diag, abovediag, belowdiag, rhs, solution);

        // boundary conditions
        gsl_vector_set(solution, 0, 0.0);
        gsl_vector_set(solution, N, 0.0);
        gsl_matrix_set_row(uallIMPL, t, solution);

        // replace old time-solution with new
        for (int i = 0; i<= N; i++) {
            gsl_vector_set( rhs, i,  gsl_vector_get( solution, i) );
        }
    }
    
    // print solution?
    for (int i=0; i<=N; i++) {
        printf("u_%d = %g\n", i,  gsl_vector_get( solution, i));
    }
    
    // // Save to file
    // 
    // Need monster-arrays with: 1) x-values, 2) computed u

    double* XVALS = new double[N+1];
    double* COMPU = new double[tsteps*(N+1)];

    for (int i=0; i<=N; i++) {
        // iterate through each x-element
        XVALS[i] = i*delta_x;
    }
    for (int t=0; t<tsteps; t++) {
        // iterate through each timestep
        for (int i=0; i<=N; i++) {
            COMPU[i + t*(N+1)] = gsl_matrix_get(uallIMPL, t, i);
        }
    }

    outputFile2(N+1, tsteps, XVALS, COMPU, &filename[0]); 
    // // // Done saving


    // HOUSEKEEPING
    gsl_vector_free(diag);
    gsl_vector_free(abovediag);
    gsl_vector_free(belowdiag);
    gsl_vector_free(rhs);
    gsl_vector_free(solution);
    gsl_matrix_free(uallIMPL);
    delete []   XVALS;
    delete []   COMPU;

}

void crank_nicolson(int N, int tsteps, double delta_x, double alpha, char* filename[]) {
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
    gsl_matrix* uallCN      = gsl_matrix_alloc(tsteps,N+1);


    // Diagonal matrix: elements
    a = c = -alpha;
    b = 2 + 2*alpha;

    for (int i=0; i < N; i++) { 
        // fill values: above and below diag
        gsl_vector_set( abovediag, i, a);
        gsl_vector_set( belowdiag, i, c);
    }
    for (int i=0; i <= N; i++) {
        // Fill diag. values
        gsl_vector_set( diag, i, b);
    }


    // Boundary conditions:
    //      u(0,t) = 1 , t>0
    //      u(d,t) = 0 , t>0

    gsl_vector_set( solution, 0, 0.0);
    gsl_vector_set( solution, N, 0.0);
    gsl_vector_set( rhs, N, 0.0);
    gsl_vector_set( rhs, 0, 1.0);
    gsl_matrix_set_row(uallCN, 0, solution);
    
    printf("abovediagN-1 = %g\n", gsl_vector_get( belowdiag, 0));

    //  Initialisation
    for (int i=1; i < N; i++) {
        // init condition
        gsl_vector_set( solution, i, initcond( i*delta_x ) );
    }
    
    // Time iteration
    for (int t = 1; t< tsteps; t++) {
        for (int i = 1; i< N; i++) {
            gsl_vector_set( rhs, i, \
                  alpha         *gsl_vector_get( solution, i-1) \
                + (2 - 2*alpha) *gsl_vector_get( solution, i) \
                + alpha         *gsl_vector_get( solution, i+1) \
                );
            }

        gsl_vector_set( rhs, 0, 0.0);       // necessary?
        gsl_vector_set( rhs, N, 0.0);       // ?

        // Solve tridiag matrix, modifies vector: solution
        gsl_linalg_solve_tridiag(diag, abovediag, belowdiag, rhs, solution);

        gsl_vector_set( solution, 0, 0.0);
        gsl_vector_set( solution, N, 0.0);
        gsl_matrix_set_row(uallCN, t, solution);
        
        }

    // print solution?
    printf("CRANK-NICOLSON SCHEME\n");
    for (int i=0; i<=N; i++) {
        printf("%g\t,", gsl_vector_get( solution, i));
    }
    printf("\n");

    // // Save to file
    // 
    // Need monster-arrays with: 1) x-values, 2) computed u

    double* XVALS = new double[N+1];
    double* COMPU = new double[tsteps*(N+1)];

    for (int i=0; i<=N; i++) {
        // iterate through each x-element
        XVALS[i] = i*delta_x;
    }
    for (int t=0; t<tsteps; t++) {
        // iterate through each timestep
        for (int i=0; i<=N; i++) {
            COMPU[i + t*(N+1)] = gsl_matrix_get(uallCN, t, i);
        }
    }

    outputFile2(N+1, tsteps, XVALS, COMPU, &filename[0]); 

    // // // Done saving


    // HOUSEKEEPING
    gsl_vector_free(diag);
    gsl_vector_free(abovediag);
    gsl_vector_free(belowdiag);
    gsl_vector_free(rhs);
    gsl_vector_free(r);
    gsl_vector_free(solution);
    gsl_matrix_free(uallCN);
    delete []   XVALS;
    delete []   COMPU;
}

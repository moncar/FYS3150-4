#include <iostream>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lib.h"
#include "headers.h"


// MHJ: Hermite zeros and weights
#define EPS 3.0e-14
#define PIM4 0.7511255444649425
#define MAXIT 10


using namespace std;


double wavefunc(double x1, double y1, double z1, \
                double x2, double y2, double z2) {
    
    double alpha = 1.0;

    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
    double r1r2=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));


    if (r1r2 > 1e-6) {
        return exp( - (alpha*alpha) * (r1*r1 + r2*r2) ) / (r1r2) ;
    }
    else {
        return 0.0;
    }
}

double wavefunc_red(double x1, double y1, double z1, \
                    double x2, double y2, double z2) {
    
    double alpha = 1.0;

    double r1r2=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

    if (r1r2 > 1e-6) {
        return 1.0 / (r1r2) ;
    }
    else {
        return 0.0;
    }
}

double wavefuncsqT1(double* x, double alpha) 
{
    
    double r1sq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r2sq = x[3]*x[3] + x[4]*x[4] + x[5]*x[5];
    //double r1r2 = sqrt((x[0]-x[3])*(x[0]-x[3]) + (x[1]-x[4])*(x[1]-x[4]) + (x[2]-x[5])*(x[2]-x[5]));

    return exp(- alpha*alpha * (r1sq + r2sq ));
}


void GaussHermite(double *x,double *w, int n)
{
    int i,its,j,m;
    double p1,p2,p3,pp,z,z1;

    m=(n+1)/2;
    for (i=1;i<=m;i++) {
        if (i == 1) {
                z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
        } else if (i == 2) {
                z -= 1.14*pow((double)n,0.426)/z;
        } else if (i == 3) {
                z=1.86*z-0.86*x[1];
        } else if (i == 4) {
                z=1.91*z-0.91*x[2];
        } else {
                z=2.0*z-x[i-2];
        }
        for (its=1;its<=MAXIT;its++) {
            p1=PIM4;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
            }
            pp=sqrt((double)2*n)*p2;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "MAXIT reached!" << endl;
        x[i]=z;
        x[n+1-i] = -z;
        w[i]=2.0/(pp*pp);
        w[n+1-i]=w[i];
    }
}

void RandNoGenGauss(double **randnos, double sigma, int dims, int N)
{
    // Calculate random numbers
    // GNU Scientific library implementation
    const gsl_rng_type * T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_taus;    // generator function
    r = gsl_rng_alloc(T);   // random number generator
    
    // C++2011 implementation:
    //random_device generator;
    //normal_distribution<double> distribution(0.0, 2.0);

    // Initialise arrays
    #pragma omp parallel for
    for (int i=0; i < dims; i++) {
        randnos[i] = new double[N];
    }
    cout << "Done initialising." << endl;

    // Calculate RVs using numbers from "generator" w/ "distribution".
    for (int j=0; j < dims; j++) {
        #pragma omp parallel for
        for (int i=0; i < N; i++) {
            randnos[j][i]    = gsl_ran_gaussian(r, sigma);
    //        x4[j][i]   = distribution(generator);
        }
        cout << "Done w/ round " << j << endl;
    }

}

double LocalEnergy1(double x1, double y1, double z1 \
                , double x2, double y2, double z2 \
                , double alpha, int N)
{
    // Calculates the local energy for the test wave function
    // PSI_T1
    
    double alfalfa = alpha*alpha;
    double r1sq;
    double r2sq;
    double r1r2; 
    
    for (int i=0; i < N; i++) {
        // loop over all steps

        r1sq = x1*x1 + y1*y1 + z1*z1;
        r2sq = x2*x2 + y2*y2 + z2*z2;
        r1r2 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

        //cout << r1sq+r2sq << " \t" << 1./r1r2 << endl;
        // Calculate local energy, store value in outputarray
        return 0.5 * alfalfa * (2. - alfalfa * (r1sq + r2sq)) \
              +0.5 * (r1sq + r2sq) + 1./r1r2;
    } 

}


int main1(int argc, char* argv[], int N) {
    
    cout << "Launching Gauss-Legendre integrator. N=" << (int)N << endl;

    // Gauss-Legendre:  integrate wave-function
    //                  which is approx zero at r=5

    double lima = -5.0;
    double limb = +5.0;


    double* x = new double[N];  // zeros/ mesh points
    double* w = new double[N];  // weights
    double intsum = 0.0;        // integration sum

    // init: mesh points and weights
    gauleg(lima, limb, x,w, N);

    // integrate
    #pragma omp parallel for reduction(+:intsum)
    for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
    for (int k=0; k < N; k++) {
    for (int l=0; l < N; l++) {
    for (int m=0; m < N; m++) {
    for (int n=0; n < N; n++) {
        intsum += w[i]*w[j]*w[k]*w[l]*w[m]*w[n] \
                * wavefunc(x[i], x[j], x[k], x[l], x[m], x[n]);
    }}}}}}

    printf("Integral sum: %g , with N=%d \n", intsum, N);
    
    // Housekeeping
    delete [] x;
    delete [] w;
    
    
    return 0;
}

int main2(int argc, char* argv[], int N) {
    
    cout << "Launching Gauss-Hermite integrator. N=" << (int)N << endl;

    // Gauss-Hermite:   integrate wave-function
    //                  using ortho func: Hermite, +/- infty
    
    double* xx = new double[N+1]; // zeros/mesh points
    double* ww = new double[N+1]; // weights
    double intsum2= 0.0;        // integration sum: G-H

    // init: mesh points and weights
    GaussHermite(xx, ww, N);

    cout << ww[0] << " " << ww[1] << endl;

    // integrate
    #pragma omp parallel for reduction(+:intsum2)
    for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
    for (int k=0; k < N; k++) {
    for (int l=0; l < N; l++) {
    for (int m=0; m < N; m++) {
    for (int n=0; n < N; n++) {
        intsum2 += ww[i]*ww[j]*ww[k]*ww[l]*ww[m]*ww[n] \
                * wavefunc_red(xx[i], xx[j], xx[k], xx[l], xx[m], xx[n]);
    }}}}}}

    printf("Integral sum: %.12g , with N=%d \n", intsum2, N);
    
    // Housekeeping
    delete [] xx;  
    delete [] ww;

    return 0;

}


int main3(int argc, char* argv[], int N) {
    
    // Monte-Carlo
    
    cout << "Launching Monte Carlo integrator. N=" << N << endl;


    // Brute-force Monte Carlo

    // init random seed 
    //srand(time(0));
    
    // Limits:
    double limit = 5.0;

    // GNU Scientific library implementation
    const gsl_rng_type * T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;    // generator function
    r = gsl_rng_alloc(T);   // random number generator

    // C++2011 implementation:
    // random_device rd;
    // uniform_real_distribution<double> uniform_real(-5.0,5.0);
   
    // Jacobi determinant
    // length[-a,a] = 2a
    double jacobiDet  = pow((2*limit), 6.0);

    double** xxx = new double*[6];
    double*  fx  = new double[N];
    #pragma omp parallel for
    for (int i=0; i < 6; i++) {
        xxx[i] = new double[N];
    }

    for (int j=0; j < 6; j++) {
        #pragma omp parallel for
        for (int i=0; i < N; i++) {
            xxx[j][i]   =  gsl_ran_flat(r, -limit, +limit); 
                        //  C++2011:
                        //= uniform_real(rd);
        }
    }


    double intsum3 = 0.0;
    double sigmasq = 0.0;
    double variance= 0.0;

    for (int j=0; j < 6; j++) {
        cout << xxx[j][0] << "  \t" << xxx[j][1] << endl;
    }

    // Calculate f(x) at the given points
    #pragma omp parallel for
    for (int i=0; i < N; i++) {
        fx[i]    = wavefunc(xxx[0][i], xxx[1][i], xxx[2][i] \
                            , xxx[3][i], xxx[4][i], xxx[5][i]);
    }

    #pragma omp parallel for reduction(+:intsum3,sigmasq)
    for (int i=0; i < N; i++) {
        intsum3 += fx[i];
        sigmasq += fx[i]*fx[i];
    }

    double integral = intsum3/N;
    sigmasq         = sigmasq/N;
    variance        = sigmasq - integral*integral;
    variance        = variance / N ; // ???

    // Scale
    integral *= jacobiDet;
    double stddev = sqrt(variance) * jacobiDet;

    printf("Integral sum: %.12g , with N=%d \n", integral, N);
    cout << "\t" << sigmasq/N << " x^2: " << integral*integral/jacobiDet << endl;
    printf("\tVariance: %.12g,\n\tstd.dev: %.12g\n", variance, stddev);

    if (argc > 3) { 
        outputFile2(N, 1, xxx[1], fx, &argv[3]);
    }
    
    // Housekeeping
    gsl_rng_free (r);
    delete [] xxx;
    delete [] fx;

    return 0;
}

int main4(int argc, char* argv[], int N) {

    // Monte Carlo: Importance sampling

    // Need: - timer for all methods
    //       - variance and stddev for all methods

    // Initialise
    double** x4 = new double*[6];   // holding: random numbers
    double*  fx4= new double[N];    // holding: function values
    clock_t t;                      // timer
    t = clock();

    // Set std.dev for Gaussian function
    double sigma = 1./sqrt(2);
    // Set Jacobi-det for integral, following MHJ
    double jacobiDet = pow(4*atan(1), 3);


//    RandNoGenGauss(x4, sigma, 6, N);

    // Calculate random numbers
    // GNU Scientific library implementation
    const gsl_rng_type * T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_taus;    // generator function
    r = gsl_rng_alloc(T);   // random number generator
    
    // C++2011 implementation:
    //random_device generator;
    //normal_distribution<double> distribution(0.0, 2.0);

    // Initialise arrays
    #pragma omp parallel for
    for (int i=0; i < 6; i++) {
        x4[i] = new double[N];
    }
    cout << "Done initialising." << endl;

    // Calculate RVs using numbers from "generator" w/ "distribution".
    for (int j=0; j < 6; j++) {
        #pragma omp parallel for
        for (int i=0; i < N; i++) {
            x4[j][i]    = gsl_ran_gaussian(r, sigma);
    //        x4[j][i]   = distribution(generator);
        }
        cout << "Done w/ round " << j << endl;
    }


    double intsum4 = 0.0;
    double sigmasq = 0.0;
    double variance= 0.0;

    for (int j=0; j < 6; j++) {
        cout << x4[j][0] << "  \t" << x4[j][1] << endl;
    }

    // Calculate f(x) at the given points
    #pragma omp parallel for
    for (int i=0; i < N; i++) {
        fx4[i]   = wavefunc_red(x4[0][i], x4[1][i], x4[2][i] \
                            , x4[3][i], x4[4][i], x4[5][i]);
    }

    cout << "Done calculating wavefunc" << endl;

    #pragma omp parallel for reduction(+:intsum4,sigmasq)
    for (int i=0; i < N; i++) {
        intsum4 += fx4[i];
        sigmasq += fx4[i]*fx4[i];
    }

    cout << "Done w/ reduction" << endl;

    double integral = intsum4/N;
    sigmasq         = sigmasq/N;
    variance        = sigmasq - integral*integral;
    variance        = variance / N ; // ???

    // Scaling
    integral        = integral * jacobiDet;
    double stddev   = sqrt(variance) * jacobiDet;

    printf("Integral sum: %.12g , with N=%d \n", integral, N);
    printf("\tStd.dev: %.12g\n", stddev);
    

    // Time spent?
    t = clock() - t;
    cout << "Elapsed time: " << (float)t/CLOCKS_PER_SEC << "s" << endl;

    // Save if argument says so
    if (argc > 3) { 
        outputFile2(N, 1, x4[0], fx4, &argv[3]);
    }

    // Housekeeping

    delete [] x4;
    delete [] fx4;

    return 0;    

}

int mainVMC1(int argc, char* argv[], int N) {
    // Variational part: using test wave function PSI_T1
    // Solution strategy:
    //  - Create 6 random nos - one for each dim
    //  - Calc wave func
    //  - Check ratio: new wavefunc sq / old wavefunc sq >= ran[0,1]?
    //  - YES:  use random nos func as starting point,
    //          save wavefunc value,
    //          save random nos as locations,
    //          iteration++
    //  - NO: keep old ones, try new sample, iteration = iteration

    // Initialise
    double** xVMC1 = new double*[N];   // holding: random numbers
    double*  fxVMC1= new double[2*N];    // holding: function values
    double*  ergLoc= new double[N];    // holding: local energy values
    clock_t t;                      // timer
    t = clock();

    // Set std.dev for Gaussian function
    double sigma = 1./sqrt(2);
    // Set Jacobi-det for integral, following MHJ
    double jacobiDet = pow(4*atan(1), 3);

    
    // Calculate random number
    // GNU Scientific library implementation
    const gsl_rng_type * T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_taus;    // generator function
    r = gsl_rng_alloc(T);   // random number generator
    
    // C++2011 implementation:
    //random_device generator;
    //normal_distribution<double> distribution(0.0, 2.0);

    // Initialise arrays
    #pragma omp parallel for
    for (int j=0; j < 6; j++) {
        xVMC1[j] = new double[N];
    }
    cout << "Done initialising." << endl;

    // Set initial trial wavefunc
    //

    //for (int j=0; j<6; j++) { xVMC1[j][0] = xi[j]; }


    
    int variations= 20;     // alpha variations
    double alphast= 0.01;    
    double alpha0 = 0.95;
    double* alphas= new double[variations];
    double* ergtot= new double[variations];
    double* stddev= new double[variations];
    double* variances = new double[variations];

    for (int a=0; a < variations; a++) { alphas[a] = alpha0 + a*alphast; }

    double  step = 0.0005;
    int     n;
    int     randcord;
    double  randno;
    double  randno2;
    double  randno3;
    double  probval;
    double  ratio;
    double  ergloc;
    double  ergloc_tmp;
    double  ergloc_sq;
    
    
    // Loop over all alphas
    #pragma omp parallel for private(ergloc,ergloc_tmp,ergloc_sq,n)
    for (int a=0; a< variations; a++) {
        cout << "Iterating for alpha=" << alphas[a] << endl;
        ergloc      = 0.0;
        ergloc_tmp  = 0.0;
        ergloc_sq   = 0.0;
        n       = 1;
        
        // Choose starting position
        double xi[6] = {1.0, 1.0, 1.0 \
                        ,2.0, 2.0, 2.0 };
        double xf[6] = {0.0, 0.0, 0.0 \
                    ,0.0, 0.0, 0.0 };
        fxVMC1[0] = wavefuncsqT1(xi, alphas[a]);


        while (n < N) {
            // Pick random coordinate [0,5] and calc ratio
            randcord = (int)floor(6*gsl_rng_uniform(r));
            // Calculate RVs
            randno = gsl_rng_uniform(r);    // in range [0,1).
            randno3= gsl_rng_uniform(r);    // in range [0,1).
            // Transform no -> [-1,1):
            randno2= 2*randno - 1.0;

            // Calculate new step:
            for (int j=0; j<6; j++) { xf[j] = xi[j]; }
            xi[randcord] = xf[randcord] + randno2*step;

            // Calc wavefunc ratio
            probval = wavefuncsqT1(xi, alphas[a]);
            ratio = probval \
                    /fxVMC1[n-1];

            // if ratio > randno
            // continue
            if (ratio > randno3) {
                for (int j=0; j<6; j++) { 
                    xVMC1[j][n] = xf[j];
                    xi[j] = xf[j];
                }
                // Calculate local energy
                ergloc_tmp= LocalEnergy1(xf[0], xf[1], xf[2] \
                                       , xf[3], xf[4], xf[5] \
                                       , alphas[a], N); 
                ergloc += ergloc_tmp;
                ergloc_sq += ergloc_tmp * ergloc_tmp;
                // Calculate probability value
                fxVMC1[n]= probval;
                n++;
            } 
        }
        // Calculate Hamiltonian 
        ergtot[a]   = 1./N * ergloc;

        // Calculate variance
        ergloc_sq   = 1./N * ergloc_sq;
        variances[a]= ergloc_sq - ergtot[a]*ergtot[a];
        stddev[a]   = sqrt(variances[a]/N);

        cout << "\tVar: " << variances[a] << ", ergloc_sq: " << ergloc_sq << endl;
        cout << "\tergloc: " << ergloc << endl;

        cout << "\tTotal energy: " << ergtot[a] << endl;
        cout << "\tStd.dev: " << stddev[a] << endl;
    }


    // Time spent?
    t = clock() - t;
    cout << "Elapsed time: " << (float)t/CLOCKS_PER_SEC << "s" << endl;
    
    double* xVMCvals = new double[N];
    for (int i=0; i<N; i++) { 
        xVMCvals[i] = xVMC1[0][i]; 
    }
//    for (int i=N; i<2*N; i++) {
//        fxVMC1[i] = ergLoc[i-N];
//    }

    if (argc > 3) { 
//        outputFile2(N, 1, xVMCvals, fxVMC1, &argv[3]);
        outputFile2(variations, 1, alphas, stddev, &argv[3]);        
    }

    // Housekeeping
    delete [] alphas;
    delete [] xVMC1;
    delete [] fxVMC1;
    delete [] xVMCvals;
    delete [] ergLoc;
    delete [] ergtot; 
    delete [] stddev;
    delete [] variances;

    return 0;    
}

    


int main(int argc, char* argv[]) {

    int N;

    if (argc < 2) {
        N = 2*2;
        cout << "Launching integrator. N=" << (int)N << endl;
    }
    else {
        N = atoi(argv[1]);
        cout << "Launching integrator. N=" << (int)N << endl;
    }

    // Display number of possible threads
    int numthreads = omp_get_max_threads();
    int numcores = omp_get_num_procs();
    cout << "Number of avail threads: " << numthreads << endl;
    cout << "Number of avail cores: " << numcores << endl;
//    omp_set_dynamic(0);
//    omp_set_num_threads(numthreads);


    cout << "[1]:\tGauss-Legendre quadrature,"  << endl;
    cout << "[2]:\tGauss-Hermite quadrature,"   << endl;
    cout << "[3]:\tBrute force Monte Carlo,"    << endl;
    cout << "[4]:\tImportance sampling Monte Carlo,"<< endl;
    cout << "[5]:\tVariational Monte Carlo, wavefunc 1,"<< endl;
    int choice;
    if (argc < 3) {
        cin >> choice;
    }
    else {
        choice = atoi(argv[2]);
    }
    switch(choice) {
        case 1: 
            cout << "Entered 1" << endl;
            int statusMain1;
            statusMain1 = main1(argc, argv, N);

            return statusMain1;
        case 2:
            cout << "entered 2" << endl;

            int statusMain2;
            statusMain2 = main2(argc, argv, N);

            return statusMain2;
        case 3:
            cout << "entered 3" << endl;

            int statusMain3;
            statusMain3 = main3(argc, argv, N);

            return statusMain3;
        case 4:
            cout << "entered 4" << endl;

            int statusMain4;
            statusMain4 = main4(argc, argv, N);

            return statusMain4;
        case 5:
            cout << "entered 5" << endl;

            int statusMain5;
            statusMain5 = mainVMC1(argc, argv, N);

            return statusMain4;
        default:
            cout << "Enter valid number..." << endl;
            return 1;
    }


}

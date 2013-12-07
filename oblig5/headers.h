#ifndef outputFile_H
#define outputFile_H
void outputFile(int N,int cols,double *fX,double *fCalcX, char *filename[]);
#endif

#ifndef outputFile2_H
#define outputFile2_H
void outputFile2(int N,int rows,double *fX,double *fCalcX, char *filename[]);
#endif

// Mains
#ifndef mainGaussLegendre_H
#define mainGaussLegendre_H
int mainGaussLegendre(int argc, char* argv[], int N);
#endif

#ifndef mainGaussHermite_H
#define mainGaussHermite_H
int mainGaussHermite(int argc, char* argv[], int N);
#endif

#ifndef mainMonteCarloBruteForce_H
#define mainMonteCarloBruteForce_H
int mainMonteCarloBruteForce(int argc, char* argv[], int N);
#endif

#ifndef mainMonteCarloImportanceSampling_H
#define mainMonteCarloImportanceSampling_H
int mainMonteCarloImportanceSampling(int argc, char* argv[], int N);
#endif

#ifndef mainMonteCarloVMC1_H
#define mainMonteCarloVMC1_H
int mainMonteCarloVMC1(int argc, char* argv[], int N);
#endif

#ifndef mainTESTSUITE_H
#define mainTESTSUITE_H
int mainTESTSUITE(int argc, char* argv[], int N);
#endif


// wavefuncs.cpp

#ifndef wavefunc_H
#define wavefunc_H
double wavefunc(double x1, double y1, double z1, \
                double x2, double y2, double z2); 

#endif

#ifndef wavefunc_red_H
#define wavefunc_red_H
double wavefunc_red(double x1, double y1, double z1, \
                    double x2, double y2, double z2);

#endif

#ifndef wavefuncT1sq_H
#define wavefuncT1sq_H
double wavefuncT1sq(double* x, double alpha, double);

#endif

#ifndef wavefuncT1_H
#define wavefuncT1_H
double wavefuncT1(double* x, double alpha, double);
#endif

#ifndef wavefuncT2sq_H
#define wavefuncT2sq_H
double wavefuncT2sq(double* x, double alpha, double beta);

#endif

#ifndef wavefuncT2_H
#define wavefuncT2_H
double wavefuncT2(double* x, double alpha, double beta);
#endif


// Misc functions

#ifndef RandNoGenGauss_H
#define RandNoGenGauss_H
void RandNoGenGauss(double **randnos, double sigma, int dims, int N);
#endif

#ifndef LocalEnergy1_H
#define LocalEnergy1_H
double LocalEnergy1(double x1, double y1, double z1 \
                , double x2, double y2, double z2 \
                , double alpha, int N);
#endif

#ifndef LocalEnergyNumerical_H
#define LocalEnergyNumerical_H
void LocalEnergyNumerical(double** x \
                        ,   double (*wvfunc)(double*, double, double) \
                        ,   double (*potential)(double*) \
                        ,   double xstep \
                        ,   int N   \
                        ,   double alpha \ 
                        ,   double beta \
                        ,   double* energies );
#endif

#ifndef potential_H
#define potential_H
double potential(double* x);
#endif

    
// External functions
#ifndef GaussHermite_H
#define GaussHermite_H
void GaussHermite(double *x,double *w, int n);
#endif



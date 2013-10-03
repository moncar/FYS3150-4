#ifndef JACOBIROTATION_H
#define JACOBIROTATION_H

void jacobiRotation(arma::mat &A, int N);

#endif

#ifndef VSTARTJACOBI_H
#define VSTARTJACOBI_H

void vStartJacobi( \
        int N, \
        double fRho_min, \
        double fRho_max, \
        double omega, \
        arma::vec* eigvals);

#endif

#ifndef VSTARTEIGVAL_H
#define VSTARTEIGVAL_H

void vStartEIGVAL( \
        int N, \
        double fRho_min, \
        double fRho_max, \
        double omega, \
        arma::vec* eigvals, \
        arma::mat* eigvec);

#endif

#ifndef VSTARTTQLI_H 
#define VSTARTTQLI_H 

void vStartTQLI(\
        int N, \
        double fRho_min, \
        double fRho_max, \
        double omega, \
        arma::vec* eigvals);

#endif

#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

void outputFile(int N,int cols,double *fX,double *fCalcX, char *filename[]);

#endif

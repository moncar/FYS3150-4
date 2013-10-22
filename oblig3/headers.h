#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

void outputFile(int N,int cols,double *fX,double *fCalcX, char *filename[]);

#endif

#ifndef asystem_H
#define asystem_H

arma::vec aSystem(double dT0, arma::vec Xii);

#endif



#ifndef fvel_H
#define fvel_H

inline double fVel(double fMass, double fRadius);

#endif


#ifndef solver_EC_H
#define solver_EC_H

void solver_EC(   arma::vec (*vFuncIn)(arma::vec, arma::vec) \
                , arma::mat *X \
                , int N \
                , int M \
                , double dStep); 

#endif

#ifndef solver_RK4_H
#define solver_RK4_H

void solver_RK4(   arma::vec (*vFuncIn)(double, arma::vec) \
                , arma::mat *X \
                , int N \
                , int M \
                , double dT0 \
                , double dStep); 

#endif

#ifndef solver_DPRI_H
#define solver_DPRI_H

void solver_DPRI(   arma::vec (*vFuncIn)(double, arma::vec) \
                , arma::mat *X \
                , int N \
                , int M \
                , double dT0 \
                , double dStep); 

#endif

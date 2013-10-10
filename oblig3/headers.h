#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

void outputFile(int N,int cols,double *fX,double *fCalcX, char *filename[]);

#endif


#ifndef solver_EC_H
#define solver_EC_H

void solver_EC(   arma::vec (*vFuncIn)(arma::vec) \
                , arma::vec *X0 \
                , arma::mat *X \
                , int N); 

#endif

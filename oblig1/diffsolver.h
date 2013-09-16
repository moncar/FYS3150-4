#ifndef fSolver_H
#define fSolver_H
void solver(double *fX, double fH, double N);

void LUmethod(int N, double *fSrcX);

void outputFile(int N, double *fX, double *fCalcX, char *filename[]);

#endif

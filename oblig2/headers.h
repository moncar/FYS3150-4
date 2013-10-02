#ifndef JACOBIROTATION_H
#define JACOBIROTATION_H

void jacobiRotation(arma::mat &A, int N);

#endif

#ifndef VSTARTJACOBI_H
#define VSTARTJACOBI_H

void vStartJacobi(double N, double fRho_min, double fRho_max);

#endif

#ifndef VSTARTEIGVAL_H
#define VSTARTEIGVAL_H

void vStartEIGVAL(double N, double fRho_min, double fRho_max);

#endif

#ifndef VSTARTTQLI_H 
#define VSTARTTQLI_H 

void vStartTQLI(int N, double fRho_min, double fRho_max);

#endif

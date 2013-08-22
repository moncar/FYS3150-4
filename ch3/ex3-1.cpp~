// Ex. 3.1
// Differentiate inverse tangent function
// Thu 13-08-22

#include <iostream>
#include <cmath>
#include <stdio.h>

using namespace std;

double fFirstDeriv(double fX, double h, double (*fFuncIn)(double) )
// Differentiation method A
// accepts function *fFuncIn, step size h and evaluation point fX
// returns value of differentiated function at fX
{
   double fRes;
   fRes = ((*fFuncIn)(fX + h) - (*fFuncIn)(fX))/h; 
   return fRes;
}

double fFirstDeriv2nd(double fX, double h, double (*fFuncIn)(double) )
// Differentiation method B
// accepts function *fFuncIn, step size h and evaluation point fX
// returns value of differentiated function at fX
{
   double fRes, fFhplus, fFhminus;
   fFhplus = (*fFuncIn)(fX + h);
   fFhminus= (*fFuncIn)(fX - h);
   fRes = (fFhplus - fFhminus) / (2*h);
   return fRes;
}

double fTestFunc(double fX)
{
   float fOut;
   fOut = (fX*fX + 3.0f)*2.0f;
   return fOut;
}

int main() 
{
   //double fRes1 = fFirst_deriv(fTestFunc(2.0));
   double fX, h;
   int N;
   // Loop: decrease step size h (and store h in array)
   //       and evaluate function at point fX
   //       using two different differentiation methods
   cout << "Enter desired number of iterations: N=";
   cin >> N;

   // array w/ dynamic size N, to hold h:
   double *fHs = new double[N];
   float  *sHs = new float[N];
   
   // array to hold results from diff. method A and B
   double *fResultsA = new double[N];
   float  *sResultsA = new float[N];
   double *fResultsB = new double[N];
   float  *sResultsB = new float[N];

   cout << endl << "Enter desired initial step size: h=";
   cin >> h;
   cout << endl << "Enter desired point to evaluate function: fX=";
   cin >> fX;

   cout << "  h  \t\t  met.A double \t  met.A single \t  met.B double \t  met.B single  " << endl;
   for (int iii=0; iii<N; iii++) {
      // calculate h
      fHs[iii] = h/(double)(iii+1);
      sHs[iii] = (float)h /(iii+1);
      // estimate derivative
      fResultsA[iii] = fFirstDeriv    (fX, fHs[iii], fTestFunc);
      sResultsA[iii] = fFirstDeriv    (fX, sHs[iii], fTestFunc);
      fResultsB[iii] = fFirstDeriv2nd (fX, fHs[iii], fTestFunc);
      sResultsB[iii] = fFirstDeriv2nd (fX, sHs[iii], fTestFunc);

      // preliminary printing to display
      printf("%e\t%e\t%e\t%e\t%e \n",fHs[iii], fResultsA[iii], \
         sResultsA[iii], fResultsB[iii], sResultsB[iii]);
   }
   
   delete [] fHs;       // housekeeping
   delete [] sHs;
   delete [] fResultsA;
   delete [] sResultsA;
   delete [] fResultsB;
   delete [] sResultsB;
   
   return 0;
}

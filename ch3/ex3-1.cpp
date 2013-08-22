// Ex. 3.1
// Differentiate inverse tangent function
// Thu 13-08-22

#include <iostream>
#include <cmath>

using namespace std;

double fFirstDeriv(double fX, double h, double (*fFuncIn)(double) )
// accepts function *fFuncIn, step size h and evaluation point fX
// returns value of differentiated function at fX
{
   double fRes;
   fRes = ((*fFuncIn)(fX + h) - (*fFuncIn)(fX))/h; 
   return fRes;
}

double fFirstDeriv2nd(double fX, double h, double (*fFuncIn)(double) )
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
   double fX, h, fRes1, fRes2;
   fX= 2.00f;
   h = 0.001;
   fRes1 = fFirstDeriv(fX, h, fTestFunc);
   fRes2 = fFirstDeriv2nd(fX, h, fTestFunc);
   cout << fRes1 << endl;
   cout << fRes2 << endl;
   return 0;
}

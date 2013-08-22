// Ex. 3.1
// Differentiate inverse tangent function
// Thu 13-08-22

#include <iostream>
#include <cmath>

using namespace std;

double fFirstDeriv(double fX, double h, double (*fFuncIn)(double) )
{
   double fRes;
   fRes = ((*fFuncIn)(fX + h) - (*fFuncIn)(fX))/h; 
   cout << fRes << endl;
   return fRes;
}

double fTestFunc(double fX)
{
   float fOut;
   fOut = (fX*fX + 3)*2.0;
   cout << fOut << endl;
   return fOut;
}

int main() 
{
   //double fRes1 = fFirst_deriv(fTestFunc(2.0));
   double h, fRes1;
   h = 0.01;
   fRes1 = fFirstDeriv(2.0, h, fTestFunc);
   //fRes1 = fFirstDeriv(fTestFunc(2.0));
   return 0;
}

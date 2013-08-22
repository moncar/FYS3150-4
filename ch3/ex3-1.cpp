// Ex. 3.1
// Differentiate inverse tangent function
// Thu 13-08-22

#include <iostream>

using namespace std;

double fFirstDeriv(double fX, double (*fFuncIn)(double) )
{
   double fRes = (*fFuncIn)(fX);
   fRes *= 3E-06;
   cout << fRes << endl;
   return fRes;
}

double fTestFunc(double fX)
{
   float fOut;
   fOut = (fX + 3)*2.0;
   cout << fOut << endl;
   return fOut;
}

int main() 
{
   //double fRes1 = fFirst_deriv(fTestFunc(2.0));
   double fRes1;
   fRes1 = fFirstDeriv(2.0, fTestFunc);
   //fRes1 = fFirstDeriv(fTestFunc(2.0));
   cout << fRes1 << endl;
   return 0;
}

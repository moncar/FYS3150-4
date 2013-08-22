// Ex. 3.1
// Differentiate inverse tangent function
// Thu 13-08-22

#include <iostream>

using namespace std;

double fFirst_deriv(double *fFuncIn(double fXvar) )
{
   return 0;
}

double fTestFunc(double fX)
{
   return pow(fX,3);
}

int main() 
{
   double fRes1 = fFirst_deriv(fTestFunc(2.0));
   return 0;
}

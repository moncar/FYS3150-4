// EXERCISE 2.1
// Wed 13-08-21

// Convert floating point number given in decimal representation to binary

#include <iostream>
#include <cmath>

using namespace std;

int nConverter(double fFloat_in)
{
   int fBinaryOut;
   int n, iii, jjj, kkk;
   // Accept fFloat_in
   // Find largest n in 2^n = fFloat_in
   n = (int)( log(fFloat_in)/log(2.) );
   cout << n << endl;

   // Create array w/ n+1 objects to hold powers of 2
   double *af2s = new double[n+1];
   for (iii=0; iii<=n; iii++)
   {
      cout << iii << endl;
      af2s[iii] = pow(2.0,iii);
   }

   // Look up given fFloat_in and compare with values in af2s[]
   // Then fill in binary array anBinary[jjj], where jjj=n -> jjj--
   jjj = n;
   int *anBinary = new int[n+1];

   while (jjj > -1)
   {
      if (fFloat_in >= af2s[jjj]) {
         // fFloat_in >= to largest power of 2
         anBinary[jjj] = 1;
         fFloat_in -= af2s[jjj];
         cout << jjj << endl;
         jjj--;
      }
      else {
         // fFloat_in < largest power of 2
         anBinary[jjj] = 0;
         cout << jjj << endl;
         jjj--;
      }
   }


   for (kkk=n;kkk >=0; kkk--) {
      cout << anBinary[kkk];
   }
   cout << endl;
   
   delete [] af2s;   // housekeeping
   delete [] anBinary;

   fBinaryOut = n;
   return fBinaryOut;
}

int main()
{
   int nTmp,nChosenValue;
   cout << "Hello. I swallow any positive integer base-10 number \nand spit out its binary counterpart. \nEnter your lucky number: ";

   cin >> nChosenValue;
   nTmp = nConverter(nChosenValue);
   return 0;
}

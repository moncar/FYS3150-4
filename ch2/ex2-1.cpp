// EXERCISE 2.1
// Wed 13-08-21

// Convert floating point number given in decimal representation to binary

#include <iostream>
#include <cmath>

using namespace std;

// if n is positive -> calculate binary integer part AND decimal part
// if n is negative -> calculate only decimal part
// afterwards:
//    Count decimals
//    Count integers
//    Normalize
//    Set exponent
//    Set sign
//    Return IEE754 compliant 32-bit result

int nConverterDecimals(double fFloatDecimals, int n)
// Converts decimals to binary decimals,
// accepts decimal number fFloatDecimals and its log2 value n.
{
   int iii, jjj;
   iii = 0;
   // create array to hold negative powers of 2
   double af2s[23];
   // change sign of n, to index array
   n = -n;
   cout << "n=" << n << endl;
   // Create array af2s of powers of 2,
   // array af2s begins at 2^(0) = af2s[0]
   for (iii=n;iii<23;iii++) {
      af2s[iii] = pow(2.,-iii);
   }

   // Set up array to hold binary decimals
   // and compare fFloatDecimals with largest af2s[i]
   int anBinary[23] = {0};

   jjj = n;
   while (jjj<23) {
      if (fFloatDecimals/af2s[jjj] > 1) {
         anBinary[jjj] = 1;
         cout << af2s[jjj] << " and anBinary = " << anBinary[jjj] << \
            "and ratio " << fFloatDecimals/af2s[jjj] << endl;
         fFloatDecimals -= af2s[jjj]; // subracts to find remainder
         jjj++;
      }
      else if (fFloatDecimals/af2s[jjj] < 1) {
         anBinary[jjj] = 0;
         cout << af2s[jjj] << " and anBinary = " << anBinary[jjj] << \
            "and ratio " << fFloatDecimals/af2s[jjj] << endl;
         jjj++;
      }
      else { // fFloatDecimals is divisible by 2
         cout << "HIT ELSE" << endl;
         anBinary[jjj] = 1;
         break;
      }
   }

   cout << ".";   // cosmetic
   for (int kkk=1; kkk<23; kkk++) {
      cout << anBinary[kkk];
   }
   cout << endl;
   
   return n; 
}

int nConverter(double fFloat_in, int n)
{
   int iii, jjj, kkk;
   // Accept fFloat_in

   if (n>0) {
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
   }

   return n;
}

int main21()
{
   int n, nTmp,nTmp2;
   double fChosenValue, fChosenValueDecimals;
   cout << "Hello. I swallow any positive base-10 number \nand spit out its binary counterpart. \nEnter your lucky number: ";

   cin >> fChosenValue;
   // Find largest n in 2^n = fChosenValue
   n = (int)( log2(fChosenValue) );
   cout << n << endl;
   // Send integer part to nConverter, and
   // decimal part to nConverterDecimals;
   fChosenValueDecimals = fmod(fChosenValue, 1);   // modulus operator fmod
   fChosenValue -= fChosenValueDecimals;
   
   nTmp = nConverter(fChosenValue, n);
   cout << "And the decimals..." << endl;
   n = (int)( log2(fChosenValueDecimals) );
   nTmp2= nConverterDecimals(fChosenValueDecimals, n);
   return 0;
}

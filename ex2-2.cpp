// Ex 2.2
// Thu 13-08-22

#include <iostream>

using namespace std;

int main()
{
   // Summing 1/n for
   //    a) n=1 to N
   //    b) n=n to 1
   // User specifies N at runtime
   
   int N;
   
   float fSumUp, fSumDown;

   cout << "Enter number of iterations N=";
   cin >> N;

   // a) Summing UP
   
   for (int n=1; n<N+1; n++) {
      fSumUp += 1./(double)n;
   }
   cout << "Sum up: s_up = " << fSumUp << endl;

   // b) Summing DOWN

   for (int n=N; n>0; n--) {
      fSumDown += 1./(double)n;
   }
   cout << "Sum down: s_down = " << fSumDown << endl;

   return 0;
}

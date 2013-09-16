// Function to store to file with name "filename".
// Takes an x-array, and a calculated f(x) array,
// and stores these as columns in a file

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
using namespace std;

ofstream ofile;


void outputFile(int N, double *fX, double *fCalcX, char *filename[])
{
    ofile.open(filename[0]);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(int i=0; i < N; i++)
    {
        ofile << setw(15) << setprecision(8) << fX[i];
        ofile << setw(15) << setprecision(8) << fCalcX[i] << endl;;
    }

    ofile.close();

}

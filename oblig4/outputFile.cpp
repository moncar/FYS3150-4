// Function to store to file with name "filename".
// Takes an x-array, and a calculated f(x) array,
// and stores these as columns in a file

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
using namespace std;

ofstream ofile;


void outputFile(int N, int cols,double *fX,double *fCalcX, char *filename[])
{
    cout << "File-printer module: printing to: " << filename[0] << endl;
    ofile.open(filename[0]);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for(int i=0; i < N; i++)
    {
        ofile << setw(18) << setprecision(8) << fX[i];
        for (int j=0; j<cols; j++) {
            ofile << setw(18) << setprecision(8) << fCalcX[i+j*N];
        }
        ofile << endl;
    }

    ofile.close();

}

void outputFile2(int N,int rows,double *fX,double *fCalcX, char *filename[])
{
    cout << "Ugly file-printer module: printing "<<rows<< " lines to: " \
         << filename[0] << endl;

    ofile.open(filename[0]);
    ofile << setiosflags(ios::showpoint | ios::uppercase);

    // Print x-axis: first row
    for (int i=0; i<N; i++) {
        ofile << setw(18) << setprecision(8) << fX[i];
    }
    ofile << endl;

    // print data: second and so-forth rows
    for (int r=0; r<rows; r++) {
        for (int i=0; i < N; i++) {
            ofile << setw(18) << setprecision(8) << fCalcX[i + r*N];
        }
        ofile << endl;
    }
    
    ofile.close();
}






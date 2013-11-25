#include <iostream>
#include <math.h>
#include <omp.h>

using namespace std;


double wavefunc(double x1, double y1, double z1, \
                double x2, double y2, double z2) {
    double alpha = 1.0;

    double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
    double r1r2=sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

    cout << r1 << " r2: " << r2 << " r1r2: " << r1r2 << endl;

    return exp( - (alpha*alpha) * (r1*r1 + r2*r2) ) / r1r2 ;
}

int main(int argc, char* argv[]) {


    // Gauss-Legendre: integrate wave-function

    int treid; 
    double sumtot = 0;
    
    #pragma omp parallel private(treid)
    {
        treid = omp_get_thread_num();
        std::cout << "Hello World! " << treid << "\n";
        sumtot += wavefunc(treid, treid+1, treid, treid+2, treid, treid+3);


    }
    
    std::cout << sumtot << std::endl;

}

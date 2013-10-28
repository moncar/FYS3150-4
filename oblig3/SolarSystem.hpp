
#ifndef SOLARSYSTEM_HPP
#define SOLARSYSTEM_HPP

class SolarSystem
{
private:
    // Initial values
    //  - stored in a resizeable vector, 
    //    that grows for every planet that is added
    std::vector<double> M ;     // Mass
    std::vector<double> vx;     // Velocity x-direction
    std::vector<double> vy;     // Velocity y-direction
    std::vector<double> x ;     // Position x
    std::vector<double> y ;     // Position y

    int nPlanets = 0;   // will change as heavenly objects are added
    int nIterations=0;  // to be defined by ConstructArray
    
    // To calculate COM and potential energy:
    double fMtot = 0.0;
    double COMx  = 0.0; // To hold temp centre of mass: x-coord
    double COMy  = 0.0; // To hold temp centre of mass: y-coord
    double tmpRadCOM_planet = 0.0;

public:
    void AddPlanet(double mass, double xpos, double ypos \
        ,double xvel, double yvel);
        // Adds planet to the class:
        //  - Increments the private vectors using initial conditions
        //  - Increments the planet counter

    int PlanetCounter() {
        // Returns the private number of planets
        return nPlanets;
    }

    arma::mat ConstructArray(int N);
        // Set up the main matrix to hold all velocities and positions
        //  -> this can be sent to SolveAll as argument 
    
    inline double Radius(double, double);
        // Calculates radius between two given points

    void TotMomentum(double&, double&, double&);
        // Return tot momentum for composite system from initial conditions
        // Arguments will be changed

    void PotEnergy(arma::mat& \
                        ,   arma::vec&  \
                        ,   arma::vec& \
                        ,   arma::vec& \
                        ,   int \
                        ,   int );

    void SolveAll(arma::mat&, double, double);
        // Calls chosen solver routine (Euler-Cromer, RK4, Dormand-Prince)
        // and solves

};

#endif

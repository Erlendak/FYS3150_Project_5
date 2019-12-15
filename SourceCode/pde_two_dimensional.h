#ifndef PDE_TWO_DIMENSIONAL_H
#define PDE_TWO_DIMENSIONAL_H
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <armadillo>


using namespace std;
using namespace arma;




int JacobiSolver(int nx, int ny, double dt, double dx,mat &A, double abstol){

    /*
    Jacobi solver for solving a given instance in two dimentions, to the system reache an
    equilibrium state.
    */

    // Algorithms operators.
    int MaxIterations = 500000; // Set a maximum of iterations.
    mat B = A;
    double D = dt/(dx*dx);

    // Iterativ solver.
    for(int k=0; k<MaxIterations; k++){

        // Backwards Euler.
        for(int i = 1; i<nx-1; i++){
            for(int j=1; j<ny-1; j++){
                A(i,j) = B(i,j) + D*(B(i+1,j) + B(i,j+1) -4*B(i,j) + B(i-1,j) + B(i,j-1));
            }
        }

        // Setup new time step, also calculating the difference between the seperath steps.
        double sum = 0.0;
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                sum += (B(i,j)-A(i,j))*(B(i,j)-A(i,j));
                B(i,j) = A(i,j);
            }
        }

        // Check if we reached equilibrium yet.
        if(sqrt(sum)<abstol){
            return k; // Returns the amount of iterations needed to reach equilibrium.
        }
    }
    cout<<"Equilibrium was not reached yet."<<endl;
    return MaxIterations; // If equilibrium is not reached, we need return the maxiterations.
}




void diffusjon_example_2dim(int n){

    /*
    Setup a simple example to showcase the JacobiSolver() function as it should.
    In this example we setup a system of 10x10 with a temperature of 0. Then one side
    of the boundary conditions is set to one, and we calculate the temperature distrubution
    to we reach a equilibrium state.
    */


    clock_t start, finish;   // Initialize time operators.

    double tol = 10e-8;      // Setup tolleranse for reaching equilibrium.
    double dx = 1.0/(n-1);   // Set delta x for the simulation.
    double dt = 0.01*dx*dx;  // set appropriate amount of time step.
    cout << dx << endl;
    cout << dt << endl;

    mat A =zeros<mat>(n,n); // Initialize matrix for simulation.

    // Set boundary condtiions.
    for(int i =0; i<n; i++){
         A(i,n-1) = 1.0;
    }

    // initialize file to save results
    ofstream gfile;
    string gfilename = "RESULTS/two_dimension.dat";
    gfile.open(gfilename);
    gfile << setiosflags(ios::showpoint | ios::uppercase);

    // Save initial condtitions.
    for( int j = 0; j<n; j++){
        for(int i=0; i<n; i++){
         gfile << setw(20) << setprecision(8) << A(j,i);
         }
    }
    gfile <<endl;

    // Simulation ;
    start = clock(); // Start time of simulation.
    int itcount = JacobiSolver(n,n,dt,dx,A,tol); // Run simulation
    finish = clock(); // Find time after simulation.

    // Give user information of the simulation.
    cout<<"Example simulation in 2 dimentions :\nEquilibrium reached after ; " << itcount<<
    " iterations. (Simulation time ; "<<( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)/60<<" Minutes)\n"  << endl;

    // Save results from simulation
    for( int j = 0; j<n; j++){
        for(int i=0; i<n; i++){
         gfile << setw(20) << setprecision(8) << A(j,i);
         }
    }
    gfile <<endl;
    gfile.close();
}




double Q_enriched(double x,double rho, double cp, double delta_t){

  /*
  In this function we find the correct amount of heat production at a given depth,
  in our simulation with an enriched radioactive soil. In this function we do not take
  radioactive decay into account.
  */

  // Heat production in the upper crust.
  if (x <= 0.16667 ){  // 0-20 Kilo meters scaled to the total 120 Kilo meters.
    return (1.44*pow(10,-6)*delta_t/(rho*cp)); // Degree Celcius
  }


  // Heat production in the lower crust.
  if (x <= 0.33334 ){ // 20-40 Kilo meters scaled to the total of 120 Kilo meters.
      return (0.35*pow(10,-6)*delta_t/(rho*cp)); // Degree Celcius
  }


  // Heat production int the mantle.
  if (x <=  1){ // 40-120 Kilometers scaled to the total of 120 Kilo meters.
    return ((0.05+0.5)*pow(10,-6)*delta_t/(rho*cp)); // Degree Celcius
  }
 throw "Syntax Error : Q_enriched() reached a serious problem. \nMay not be scaled proparly";
}


void simulation_with_enrichment(int nx,int ny, double dt, double dxy, mat &A, int tsteps){

    /*
    Simulate the temperature distrubution in the lithosphere to an equilibrium state, which we can use as an initial condition
    in our main simulation. We take into account the heat production in the lithosphere and a radioactive heat contribution in
    the mantle, but we do not yet take into account the decay in radio activity, since we are looking for an equilibrium state.
    */

    clock_t start, finish; // initialize time operators.


    //Simulation spesific constants ;

    const double rho = pow(3.510,3);      // Densety              ;       Kg/m^3
    const double cp = 1000;               // Heat capasity        ;       J/Kg/ degree Celcius
    const double k = 2.5;                 // Thermal conductivity ;       W/m/ degree Celcius

    //Simulation spesific operators ;
    double scaled_x;                      // Scaled lenght of x from 120 Kilo meters to 1.
    double Q;                             // The specific heat production at a given depth, time and time step.
    double abstol = 10e-4;                // The tolerance set to consider a state reached eqilibrium.

    // Initialize file to save results.
    ofstream gfile;
    string gfilename = "RESULTS/simulation_with_enrichment.dat";
    gfile.open(gfilename);
    gfile << setiosflags(ios::showpoint | ios::uppercase);

    // Save initial condition for this simulation.
    for( int i = 0; i<nx; i++){
        for(int j=0; j<ny; j++){
         gfile << setw(20) << setprecision(8) << A(i,j);
         }
    }
    gfile<<endl;

    //Backward Euler ;
    mat B = A; //Gammel A
    double D = dt/(dxy*dxy)*k/(rho*cp);

    start = clock(); // Start taking time of the simulation.
    //Iterativ solver.
    for(int k=0; k<tsteps; k++){

        for(int i = 1; i<nx-1; i++){
            // Precalculate the heat production at a given depth, time and time interval.
            scaled_x =(double)i/ (double)nx;
            Q = Q_enriched(scaled_x,rho,cp,dt) ;
            for(int j=1; j<ny-1; j++){
                A(i,j) = ( (B(i,j) + D*(B(i+1,j) + B(i,j+1) -4*B(i,j) + B(i-1,j) + B(i,j-1))))+Q;
            }
        }

        // Setup new time step, also calculating the difference between the seperath steps.
        double sum = 0.0;
        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                sum += (B(i,j)-A(i,j))*(B(i,j)-A(i,j));
                B(i,j) = A(i,j);
            }
        }
        // Check if we reached equilibrium yet.
        if(sqrt(sum)<abstol){
            // Gives information on how the simulation have behaved.
            finish = clock(); // Gives the time at which the simulation has found equilibrium state.
            cout<<"Simulation with radioactive heat production :\nEquilibrium reached at "<<
            k*dt/60/60/24/365/pow(10,9)<<" Giga years. (Simulation time ; "<<
            ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)/60<<" Minutes)\n" <<endl;
            break; // Equilibrium reached breaks time loop and countiues to next step.
            k=tsteps; // alternativ for break, Set the iterator out of range so that this wil be last iteration.
        }
    }

    // Saves the final result when equilibrium has been reached
    for( int i = 0; i<nx; i++){
        for(int j=0; j<ny; j++){
         gfile << setw(20) << setprecision(8) << A(i,j);
         }
    }
    gfile <<endl;
    gfile.close();
}




double Q_enriched_decay(double x,double rho, double cp, double delta_t, double radio_active_contribution){

  /*
  In this function we find the correct amount of heat production at a given depth,
  in our simulation with an enriched radioactive soil, also taking into account
  the decay of radioactive substances.
  */

  // Heat production in the upper crust.
  if (x <= 0.16667){ // 0-20 Kilo meters scaled to the total 120 Kilo meters.
    return (1.44*pow(10,-6)*delta_t/(rho*cp)); // Degree Celcius
  }


 // Heat production in the lower crust.
  if (x <= 0.33334 ){ // 20-40 Kilo meters scaled to the total 120 Kilo meters.
      return (0.35*pow(10,-6)*delta_t/(rho*cp));// Degree Celcius
  }


 // Heat production int the mantle.
  if (x <= 1 ){ // 40-120 Kilo meters scaled to the total 120 Kilo meters.
    return ((0.05+(radio_active_contribution))*pow(10,-6)*delta_t/(rho*cp));// Degree Celcius
  }
  throw "Syntax Error : Q_enriched_decay() reached a serious problem. \nMay not be scaled proparly";
}


void simulation_implemented_enrichment_decay(int nx,int ny, double dt, double dxy, mat &A, int tsteps,int start_step){

    /*
    The main simulation of the project. Simulating the temperature distrubution in the lithosphere,
    with heat production and a radioactive heat contribution that is varing with time due to the radioactive substance
    is decaying at aspecific half life interval.
    */

    //Simulation spesific constants ;

    double rho = pow(3.510,3);      // Densety              ;       Kg/m^3
    double cp = 1000;               // Heat capasity        ;       J/Kg/ degree Celcius
    double k = 2.5;                 // Thermal conductivity ;       W/m/ degree Celcius

    //Simulation spesific operators
    double time;      // The spesefic simulated time.
    double U;         // radioactive heat contribution from Uran.
    double Th;        // radioactive heat contribution from Thorium, the Norwegian national element.
    double K;         // radioactive heat contribution from unstable Kalium
    double scaled_x;  // Scaled lenght of x from 120 Kilo meters to 1.
    double Q;         // The specific heat production at a given depth, time and time step.


    //Backward Euler ;
    mat B = A; //Gammel A
    double D = dt/(dxy*dxy)*k/(rho*cp);

    //iterativ solver.
    for(int k=start_step; k<tsteps; k++){

        // Precalculate radioactive contribution in the manlte.
        time = k*dt/60/60/24/365/pow(10,9);   // Time in Giga years
        U  = 0.2 *pow(0.5 , time/4.47);       //  Half life  ;  4.47  Giga years;
        Th = 0.2 *pow(0.5 , time/14.0);       //  Half life  ; 14.0   Giga years;
        K  = 0.1 *pow(0.5 , time/1.25);       //  Half life  ;  1.25  Giga years;

        for(int i = 1; i<nx-1; i++){
            // Precalculate the heat production at a given depth, time and time interval.
            scaled_x =(double)i/ (double)nx;
            Q = Q_enriched_decay(scaled_x,rho,cp,dt,(U+Th+K));
            for(int j=1; j<ny-1; j++){

                A(i,j) = ( (B(i,j) + D*(B(i+1,j) + B(i,j+1) -4*B(i,j) + B(i-1,j) + B(i,j-1))))+Q;
            }
        }

        B=A;
    }
}




#endif // PDE_TWO_DIMENSIONAL_H

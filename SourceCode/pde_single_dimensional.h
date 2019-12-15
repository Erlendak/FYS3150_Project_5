#ifndef PDE_SINGLE_DIMENSIONAL_H
#define PDE_SINGLE_DIMENSIONAL_H
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <armadillo>


using namespace std;
using namespace arma;




mat forward_euler(double delta_x, double delta_t, double lenght, double time, vec u){

    /*
    The explicit scheme using the Forward Euler method. This function takes in the lenght step size, time step size,
    the total lenght, the total lenght and initial condition. Keeping the boundary conditiion we find the distribution
    after the given time. We are returning the results as a matrix that evolves based of time.
    */

    int n = (int)(lenght/delta_x);    // Find amount of points needed for a lenght with a specific lenght step.
    int tsteps = (int)(time/delta_t); // Find amount of point needed for a time with a specific time step.

    //Setup the result matrix to keep eye on the results evloution based on time.
    mat ans(tsteps+1,n+1);
    // Saves the initial conditions.
    for (int i = 0; i <= n; i++){
    ans(0,i) = u(i);
    }

    //Forward Euler.
    vec unew = u;
    double alpha = delta_t/(delta_x*delta_x);

   // Time integration
   for (int t = 1; t <= tsteps; t++) {

     for (int i = 1; i < n; i++) {
       unew(i) = alpha * u(i-1) +( (1 - 2*alpha) * u(i) )+ (alpha * u(i+1));
     }

     // Save results at given time iteration.
     for (int i = 0; i <= n; i++){
       ans(t,i) = unew(i);
     }
     u = unew; // Set the old time step to new
   }
   return ans;
}




void fast_special_thomas_algorithm(double a,double d,double c,vec &y,vec &u, int N){

    /*
    This is the special Thomas algorithm from project 1, which solves a equation system on a tridiagonal matrix form.
    Asumming that the upper and lower diagonal is constant and that the main diagonal is constant through the equation.
    */

    double ac0 = a; // constant for upperdiagonal and lower diagonal

    vec b =  u; // Makes copy of old time step so that we could do row operations on it.
    vec tmp_c(N-2); // Saves the Gauss-row operations on the upper diagonal.

    // Initial conditions.
    b(0) = 0; // Set up boundary conditions for the equation system.
    tmp_c(0) = tmp_c(0)/d; // First Guass operation on the upper diagonal
    b(1) = b(1)/d; // Frist Gauss operation on the solution.

    //Forward substitution.
    for (int i =2; i<(N-1);i++ ){
      tmp_c(i-1) = (     ac0  /  ( d - (ac0*tmp_c(i-2))  )     );
      b(i) = ( (b(i) -(ac0*b(i-1) ) )   / (d- (ac0*tmp_c(i-2)))   );
    }

    //Backward substitution.
    for (int i = 1; i<N-1;i++){
      b((N-1)-i) = b((N-1)-i) - (tmp_c((N-2)-i)*b(N-i) );
    }
    y = b;
}


mat backward_euler(double delta_x, double delta_t, double lenght, double time, vec u){
  /*
  Implicit scheme using the backward Euler method,
  */

   int n = (int)(lenght/delta_x)+1;    // Find amount of points needed for a lenght with a specific lenght step.
   int tsteps = (int)(time/delta_t); // Find amount of point needed for a time with a specific time step.

   // Save results.
   mat ans(tsteps+1,n);
   for (int i = 0; i < n; i++){
    ans(0,i) = u(i);
   }


   // Initial conditions
   vec v(n+1); // This is u  of Av = y
   for (int i = 1; i <= n; i++) {
      v(i) = u(i-1);
   }
   vec y = v; // Right side of matrix equation Au=y, the solution at a previous step


   // Backwards Euler :
   double alpha =  delta_t/(delta_x*delta_x);
   double a = - alpha;
   double c = - alpha;
   double b = 1 + (2*alpha);

   // Time iteration
   for (int t = 1; t <= tsteps; t++) {
      //  here we solve the tridiagonal linear set of equations,
      fast_special_thomas_algorithm(a, b, c, y, v, n+1);
      // boundary conditions
      y(0) = u(0);
      y(n) = u(n-1);

      for (int i = 1; i <= n; i++){
        ans(t,(i-1)) = y(i);
      }

      // replace previous time solution with new
      for (int i = 0; i <= n; i++) {
      v(i) = y(i);
      }
      //  You may consider printing the solution at regular time intervals
      //....   // print statements
   }  // end time iteration
   //...
  //cout<<ans<<endl;
  return(ans);


}




mat crank_nicolson(int n, int tsteps, double delta_x,double alpha){
    /*
    Implicit scheme using the Crank Nicolson method,
    */
    mat ans(tsteps,n);
    double a,b,c;
    vec u(n+1);
    vec r(n+1);

    a = c = -alpha;
    b = 2 + 2*alpha;
    u(0)=1;
    for(int i=0; i<n; i++){u(i)=0; }
    for(int i=1; i<n; i++){
        ans(0,i)= alpha*u(i-1) + (2-2*alpha)*u(i) + alpha*u(i+1);
    }
    for( int t = 1; t< tsteps; t++){
        for(int i=1; i<n; i++){
            r(i) = alpha*u(i-1) + (2-2*alpha)*u(i) + alpha*u(i+1);
        }
        r(0) = 0;
        r(1) = 0;
        r(n) = 1;

        fast_special_thomas_algorithm(a,b,c,u,r,n+1);
        u(0) = 0;
        u(n) = 1;
        for(int i=1; i<=n; i++){
            ans(t,i-1)= u(i);
        }
    }
    //cout << u<<endl;
    return(ans);

}




double Q(double x,double rho, double cp, double delta_t){

  /*
  In this function we find the correct amount of heat production at a given depth,
  where we do not have a radioactive heat contribution in the mantle. So this function
  gives the appropriate initial conditions for simulating with radioactive contribution.
  */


  // Heat production in the upper crust.
  if (x <= 0.16667 ){ // 0-20 Kilo meters scaled to the total 120 Kilo meters.
    return (1.44*pow(10,-6)*delta_t/(rho*cp)); // Degree Celcius
  }


  // Heat production in the lower crust.
  if (x <= 0.33334 ){ // 20-40 Kilo meters scaled to the total 120 Kilo meters.
      return (0.35*pow(10,-6)*delta_t/(rho*cp)); // Degree Celcius
  }


 // Heat production int the mantle.
  if (x <= 1 ){ // 40-120 Kilo meters scaled to the total 120 Kilo meters.
    return (0.05*pow(10,-6)*delta_t/(rho*cp)); // Degree Celcius
  }
  throw "Syntax Error : Q() reached a serious problem. \nMay not be scaled proparly";
}


vec simulation_no_radioactive_enrichment(int n, int tsteps, double delta_x, double delta_t){

    /*
    The first simulation of the physical system, we are simulating the heat distrubution in the lithosphere.
    Where the heat production is constant at a specific depth so we are expecting some sort of symmetry, so
    we are only simulating in one dimention. The result from this simulation can we use further on as
    bounduary conditions, for the simulation with a specific radioactive heat contribution in a limited area.
    */

     clock_t start, finish; // initialize time operators.


    //Simulation spesific constants ;

    double temp_mantle = 1300;      // Mantle temperature   ;  1300 degree Celsius.
    double temp_crust  = 8;         // Surface temperature  ;     8 degree celsius.
    double rho = pow(3.510,3);      // Densety              ;       Kg/m^3
    double cp = 1000;               // Heat capasity        ;       J/Kg/ degree Celcius
    double k = 2.5;                 // Thermal conductivity ;       W/m/ degree Celcius

    //Simulation spesific operators ;
    double scale_x;                 // Scaled lenght of x from 120 Kilo meters to 1.
    double sum;                     // Finds the total diffrence between the diffrent time steps.
    double tol = 10e-8;             // The tolerance set to consider a state reached eqilibrium.


    //Backward Euler ;
    double alpha = delta_t/(delta_x*delta_x); // Backwards Euler specific constant.
    double a = - alpha*k/(rho*cp); // Backwards Euler specific operators.
    double c = - alpha*k/(rho*cp); // Backwards Euler specific operators.
    double b = 1 + (2*alpha)*k/(rho*cp); // Backwards Euler specific operators.
    vec u(n+1); // This is u  of Au = y
    vec y(n+1); // Right side of matrix equation Au=y, the solution at a previous step

    // Initial conditions
    for (int i = 1; i < n; i++) {
      y(i) = u(i) = i*temp_mantle/n; // Initial conditions, set so that we should reach equalibrium state the fastest.
    }

    // Boundary conditions
    y(n) = u(n) = temp_mantle;
    u(0) = y(0)= temp_crust;


    start = clock(); // Start taking time of the simulation.
    // Start time iteration.
    for (int t = 1; t <= tsteps; t++) {

      // Solve the tridiagonal linear set of equations using the fast special Thomas algorithm.
      fast_special_thomas_algorithm(a, b, c, y, u, n+1);

      // Take heat production into account at.
      for (int i = 1; i < n; i++) {
        scale_x = (double)i/ (double)n;
        y(i) = y(i)+(Q(scale_x,rho,cp, delta_t));
      }

      //Check up upon that the boundary conditions is set corrected.
      u(0)=y(0) = temp_crust;
      u(n)=y(n) = temp_mantle;

      // Setup new time step, also calculating the difference between the seperath steps.
      sum = 0.0;
      for (int i = 0; i <= n; i++) {
        sum += (u(i)-y(i)) * ( y(i)- u(i) );
        u(i) = y(i); // Replace previous time solution with new.
      }

      // Check if we reached equilibrium yet.
      if(sqrt(sum)<tol){
        // Gives information on how the simulation have behaved.
        finish = clock(); // Gives the time at which the simulation has found equilibrium state.
        cout<<"Simulation in one dimention:\nEquilibrium reached at "<<t*delta_t/60/60/24/365/pow(10,9)<<
        " Giga years. (Simulation time ; "<< ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)/60<<" minutes)\n" <<endl;
        return y; // Returns the eqilibrium state of the system. for further use.
      }

    }// The end of time iteration.

    throw "Simulation Error : simulation_no_radioactive_enrichment() \nNever reached the eqilibrium state.";
}




#endif // PDE_SINGLE_DIMENSIONAL_H

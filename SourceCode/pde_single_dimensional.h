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

inline double g(double x){return 0;} //abs(0.01*(x-0.5));}

inline double a(double t){return (1); //+(t/10));
}
inline double b(double t){return (0);//(2*t/20);
}

inline double f(double x){return 100.0*exp(-10.0*x);
}



void tridag(double a,double d,double c,vec &y,vec &u, int N){
/*
  double ac0 = a;
  vec b = u;

  vec tmp_c(N-2);//(N-2)

   //b(0) = 0;
   tmp_c(0) = tmp_c(0)/d;
   b(1) = b(1)/d;
   for (int i =2; i<(N-1);i++ ){
       tmp_c(i-1) = (     ac0  /  ( d - (ac0*tmp_c(i-2))  )     );
       b(i) = ( (b(i) -(ac0*b(i-1) ) )   / (d- (ac0*tmp_c(i-2)))   );// #d[i] - (_b[i-1] * a[i-1] )  );
    }

   for (int i = 1; i<N-1;i++){
       b((N-1)-i) = b((N-1)-i) - (tmp_c((N-2)-i)*b(N-i) );
   }
   //for (int i =1; i<N;i++){
    // b(i) =u(i-1);
   //}
   y=b;
*/
    double ac0 = a;
    vec b =  u;

    vec tmp_c(N-2);

     b(0) = 0;
     tmp_c(0) = tmp_c(0)/d;
     b(1) = b(1)/d;
     for (int i =2; i<(N-1);i++ ){
         tmp_c(i-1) = (     ac0  /  ( d - (ac0*tmp_c(i-2))  )     );
         b(i) = ( (b(i) -(ac0*b(i-1) ) )   / (d- (ac0*tmp_c(i-2)))   );// #d[i] - (_b[i-1] * a[i-1] )  );
      }

     for (int i = 1; i<N-1;i++){
         b((N-1)-i) = b((N-1)-i) - (tmp_c((N-2)-i)*b(N-i) );
     }
     y = b;

  }





mat forward_euler(double delta_x, double delta_t){
    int n = (int)(1/delta_x)-1;
    int tsteps = (int)(1/delta_t)-1; // delta_t/delta_x^2  =<  1/2
    mat ans(tsteps+1,n+1);

    vec u(n+1); // Setup our inital condition vector
    vec unew(n+1);

    double x;
    double alpha = delta_t/(delta_x*delta_x);
    //  First we set initialise the new and old vectors
    //  Here we have chosen the boundary conditions to be zero.
    //  n+1 is the number of mesh points in x
    //  Armadillo notation for vectors
    u(0) = unew(0) = a(0);
    u(n) = unew(n) = b(tsteps);
    for (int i = 1; i < n; i++) {
        x =  i*delta_x;
        //  initial condition
        u(i) =  g(x);
        //  intitialise the new vector
        unew(i) = 0;
    }
   for (int i = 0; i <= n; i++){
    ans(0,i) = u(i);
   }
   // Time integration
   for (int t = 1; t <= tsteps; t++) {
     unew(0) = a(t*delta_t);
     unew(n) = b(t*delta_t);
     for (int i = 1; i < n; i++) {
             // Discretized diff eq
             unew(i) = alpha * u(i-1) +( (1 - 2*alpha) * u(i) )+ (alpha * u(i+1));
      }
     for (int i = 0; i <= n; i++){
       ans(t,i) = unew(i);
     }
     u = unew;

      //cout<<unew <<endl;
       //  note that the boundaries are not changed.
     }
   //cout<< ans <<endl;
   return ans;
   //cout<< ans <<endl;
};

//  parts of the function for backward Euler

mat backward_euler(int n, int tsteps, double delta_x, double alpha)
{

   mat ans(tsteps+1,n);

   double a, b, c;
   vec u(n+1); // This is u  of Au = y
   vec y(n+1); // Right side of matrix equation Au=y, the solution at a previous step


   // Initial conditions
   for (int i = 1; i < n; i++) {
      y(i) = u(i) = g(delta_x*i);
   }
   // Boundary conditions (zero here)
   y(n) = u(n) = 1;
   u(0) = y(0)= 0;
   // Matrix A, only constants
   a = c = - alpha;
   b = 1 + (2*alpha);
   for (int i = 1; i <= n; i++){
    ans(0,(i-1)) = u(i);
   }

   // Time iteration
   for (int t = 1; t <= tsteps; t++) {
      //  here we solve the tridiagonal linear set of equations,
      tridag(a, b, c, y, u, n+1);
      // boundary conditions
      y(0) = 0;
      y(n) = 1;

      for (int i = 1; i <= n; i++){
        ans(t,(i-1)) = y(i);
      }

      // replace previous time solution with new
      for (int i = 0; i <= n; i++) {
      u(i) = y(i);
      }
      //  You may consider printing the solution at regular time intervals
      //....   // print statements
   }  // end time iteration
   //...
  //cout<<ans<<endl;
  return(ans);


}

void tridiag2(double a, double b, double c,vec r, vec &u, int n){

    vec g(n+1);
    vec d(n+1);

    for(int i = 1; i<=n; i++){
        g(i) = r(i);
    }

     d(0) = d(n) = 1;

    //Forward sub
    for(int i=1;i<=n-1;i++){
      d(i) = b-(a*c)/d(i-1);
      g(i) = g(i)-(a*g(i-1))/d(i-1);
    }

    //Backward sub
    for(int i=n-1; i>0; i--){
        u(i) = (g(i)-c*u(i+1))/d(i);
    }

}
mat crank_nicolson(int n, int tsteps, double delta_x,double alpha){
    mat ans(tsteps,n+1);
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
        r(n) = 1;

        tridiag2(a,b,c,r,u,n);
        u(0) = 0;
        u(n) = 1;
        for(int i=0; i<=n; i++){
            ans(t,i)= u(i);
        }
        }
    return(ans);
    //u.print();
}


double Q_sim1(double x,double rho, double cp, double delta_t){
  //return 0;
  double L = 120;
  //Upper crust
  if (x <= ( (1/L) *20 ) ){ // Scaled to L
    return (1.44*pow(10,-6)*delta_t/(rho*cp)); // Need to be scaled appropriatly
  }


  //Lower crust
  if (x <= ( (1/L) *40 ) ){ // Scaled to L
      return (0.35*pow(10,-6)*delta_t/(rho*cp)); // Need to be scaled appropriatly
  }


  //The mantle
  if (x <= ( (1/L) *120 ) ){ // Scaled to L
    return (0.05*pow(10,-6)*delta_t/(rho*cp)); // Need to be scaled appropriatly
  }
  cout<<"Check"<<endl;
 return 0;
}

vec simulation_before_radioactive_enrichment(int n, int tsteps, double delta_x, double delta_t){
     clock_t start, finish;


    //Simulation spesific constants ;

    double temp_mantle = 1300;      // Mantle temperature   ;  1300 degree Celsius.
    double temp_crust  = 8;         // Surface temperature  ;     8 degree celsius.
    double rho = pow(3.510,3);      // Densety              ;       Kg/m^3
    double cp = 1000;               // Heat capasity        ;       J/Kg/ degree Celcius
    double k = 2.5;                 // Thermal conductivity ;       W/m/ degree Celcius
    //Backward Euler ;

       //mat ans(tsteps+1,n);
       double scale_op, sum;
       double tol = 0.1;
       double a, b, c;
       vec u(n+1); // This is u  of Au = y
       vec y(n+1); // Right side of matrix equation Au=y, the solution at a previous step

       // Initial conditions
       for (int i = 1; i < n; i++) {
          y(i) = u(i) = i*temp_mantle/n; // Initial condtions, set for fastest equalibrium state
       }

       double alpha = delta_t/(delta_x*delta_x);
       // Boundary conditions
       y(n) = u(n) = temp_mantle;
       u(0) = y(0)= temp_crust;
       // Matrix A, only constants
       a = c = - alpha*k/(rho*cp);
       b = 1 + (2*alpha)*k/(rho*cp);

       // Save our results for time the first time step.
      // for (int i = 1; i <= n; i++){
        //ans(0,(i-1)) = u(i);
       //}
       start = clock();
       // Time iteration
       for (int t = 1; t <= tsteps; t++) {
          //  here we solve the tridiagonal linear set of equations,
          tridag(a, b, c, y, u, n+1);



          // replace previous time solution with new
          for (int i = 1; i < n; i++) {
            scale_op = (double)i/ (double)n;
            y(i) = y(i)+(Q_sim1(scale_op,rho,cp, delta_t));
           // cout<<scale_op<<endl;

          }
          u(0)=y(0) = temp_crust;
          u(n)=y(n) = temp_mantle;

          //for (int i = 1; i <= n; i++){
            //ans(t,(i-1)) = u(i);
          //}


          sum = 0.0;

          // replace previous time solution with new
          for (int i = 0; i <= n; i++) {
            sum += (u(i)-y(i)) * ( y(i)- u(i) );
            u(i) = y(i);
          }

          if(sqrt(sum)<tol){
              finish = clock();

              cout<<"Equilibrium reached at "<<t*delta_t/60/60/24/365/pow(10,9) << " Giga years. (Simulation time ; "<< ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)/60<<" minutes)" <<endl;
              return y;
          }

          //  You may consider printing the solution at regular time intervals
          //....   // print statements
       }  // end time iteration
       //...
      //cout<<ans<<endl;
      return(y);



}

#endif // PDE_SINGLE_DIMENSIONAL_H

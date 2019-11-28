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

void tridag(double a,double b,double c,vec y,vec u, int N){
    int n = N-1;
    vec tmp(N);

    u(0)= u(N-1) = 0.0;
    tmp(1) = tmp(N-1) = b;

    //Forward sub
    for(int i=2;i<=n;i++){
      tmp(i) = b-(a*c)/tmp(i-1);
      u(i) = u(i)-(a*u(i-1))/tmp(i-1);
    }

    //Backward sub
    //u[n-1] = g[n-1]/d[n-1];
    y(n-1) = u(n-1)/tmp(n-1);
    for(int i=n-2; i>0; i--){
        y(i) = (u(i)-c*y(i+1))/tmp(i);
    }


}

void thomas_algorithm(int n, double delta_x, double upper_diag, double main_diag, double lower_diag, double *u){


    double h = delta_x;

    double *x = new double[n+1]; // Position vector removebal?



    double *g = new double[n+1]; // Temperarly solution vector

    double *d = new double[n+1]; // New temporary upperdiagonal


    double *a = new double[n+1]; // Upper diagonal
    double *b = new double[n+1]; // Main diagonal
    double *c = new double[n+1]; // Lower diagonal

    for(int i = 1; i<=n; i++){
        x[i] = i*h;
        g[i] = (h*h)*f(i*h);// Need to update this one

        a[i] = upper_diag;
        b[i] = main_diag;
        c[i] = lower_diag;
    }

    u[0]= u[n] = 0.0; d[1] = d[n] = b[1];

    //Forward sub
    for(int i=2;i<=n;i++){
        d[i] = b[i]-(a[i]*c[i-1])/d[i-1];
      g[i] = g[i]-(a[i]*g[i-1])/d[i-1];
    }

    //Backward sub
    //u[n-1] = g[n-1]/d[n-1];
    u[n-1] = g[n-1]/d[n-1];
    for(int i=n-2; i>0; i--){
        u[i] = (g[i]-c[i+1]*u[i+1])/d[i];
    }

    delete[] x; delete[] g; delete[]  d; delete[] a; delete[] b; delete[] c;
  //return (u);
    }


void forward_euler(double delta_x){

    int n = 1/delta_x;

    double *u = new double[n+1]; // Final solution vector
    double delta_t;


    double upper_diag  = -1;
    double main_diag   =  2;
    double lower_diag  = -1;

   thomas_algorithm(n, delta_x, upper_diag, main_diag, lower_diag, u);
   // ofstream ofile;
    //      ofile << "       x:             approx:          exact:       relative error" << endl;
    for (int i = 1; i < n;i++) {
        cout<<u[i] <<endl;
    //    double epsilon = log10(abs((u_val(x[i])-u[i])/u_val(x[i])));
     //   ofile << setw(15) << setprecision(8) << u_val(x[i]);
      //  ofile << setw(15) << setprecision(8) << u[i];
       // ofile << setw(15) << setprecision(8) << x[i];
       // ofile << setw(15) << setprecision(8) << epsilon <<endl;
    }

   delete[] u;
};



void backward_euler(double delta_x, double delta_t){

   int n = 1/delta_x;

   double *u = new double[n+1]; // Final result vector

   double alpha = delta_t/(delta_x*delta_x);


   double upper_diag  = -alpha;
   double main_diag   =  1+(2*alpha);
   double lower_diag  = -alpha;

   thomas_algorithm(n, delta_x, upper_diag, main_diag, lower_diag, u);
   // ofstream ofile;
    //      ofile << "       x:             approx:          exact:       relative error" << endl;
    for (int i = 1; i < n;i++) {
        cout<<u[i] <<endl;
    //    double epsilon = log10(abs((u_val(x[i])-u[i])/u_val(x[i])));
     //   ofile << setw(15) << setprecision(8) << u_val(x[i]);
      //  ofile << setw(15) << setprecision(8) << u[i];
       // ofile << setw(15) << setprecision(8) << x[i];
       // ofile << setw(15) << setprecision(8) << epsilon <<endl;
    }

   delete[] u;

};

void crank_nicolson(double delta_x, double delta_t){


    int n = 1/delta_x;

    double *u = new double[n+1]; // Final solution vector

    double alpha = delta_t/(delta_x*delta_x);


    double upper_diag  = -alpha;
    double main_diag   =  1+(2*alpha);
    double lower_diag  = -alpha;

    thomas_algorithm(n, delta_x, upper_diag, main_diag, lower_diag, u);
    // ofstream ofile;
     //      ofile << "       x:             approx:          exact:       relative error" << endl;
     for (int i = 1; i < n;i++) {
         cout<<u[i] <<endl;
     //    double epsilon = log10(abs((u_val(x[i])-u[i])/u_val(x[i])));
      //   ofile << setw(15) << setprecision(8) << u_val(x[i]);
       //  ofile << setw(15) << setprecision(8) << u[i];
        // ofile << setw(15) << setprecision(8) << x[i];
        // ofile << setw(15) << setprecision(8) << epsilon <<endl;
     }

    delete[] u;
};

void barbarian_forward_euler(double delta_x, double delta_t){


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
   cout<< ans <<endl;
   //cout<< ans <<endl;
};

//  parts of the function for backward Euler
void clone_backward_euler(int n, int tsteps, double delta_x, double alpha)
{
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
   b = 1 + 2*alpha;
   // Time iteration
   for (int t = 1; t <= tsteps; t++) {
      //  here we solve the tridiagonal linear set of equations,
      tridag(a, b, c, y, u, n+1);
      // boundary conditions
      u(0) = 0;
      u(n) = 0;
      // replace previous time solution with new
      for (int i = 0; i <= n; i++) {
     y(i) = u(i);
      }
      //  You may consider printing the solution at regular time intervals
      //....   // print statements
   }  // end time iteration
   //...
}

#endif // PDE_SINGLE_DIMENSIONAL_H

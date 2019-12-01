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


/*
void crank_nicolson(int n, int tsteps, double delta_x, double alpha)
{

   mat ans(tsteps+1,n+1);

   double a, b, c;
   vec u(n+1); // This is u  of Au = y
   vec y(n+1); // Right side of matrix equation Au=y, the solution at a previous step


   // Initial conditions
   for (int i = 1; i < n; i++) {
      y(i) = u(i) = g(delta_x*i);
   }
   // Boundary conditions (zero here)
   y(n) = u(n) = 0;
   u(0) = y(0)= 1;
   // Matrix A, only constants
   a = c = - alpha;
   b = 2 + 2*alpha;
   // Time iteration
   for (int t = 1; t <= tsteps; t++) {
      //  here we solve the tridiagonal linear set of equations,
      tridag(a, b, c, y, u, n+1);
      // boundary conditions
      y(0) = 0;
      y(n) = 1;
            // replace previous time solution with new
      for (int i = 0; i <= n; i++) {
     u(i) = y(i);
      }
      for (int i = 0; i <= n; i++){
        ans(t,i) = y(i);
      }
      //  You may consider printing the solution at regular time intervals
      //....   // print statements
   }  // end time iteration
   //...
  cout<<ans<<endl;

}
*/



/*

  */
#endif // PDE_SINGLE_DIMENSIONAL_H

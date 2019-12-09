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

int JacobiSolver(int n, double dt, double dx,mat &A, double abstol){
    int MaxIterations = 50000;
    mat B = zeros<mat>(n,n); //Gammel A

    double D = dt/(dx*dx);

    //Initialbetingelse
    for(int i=1; i <n-1; i++){
        for(int j=1; j<n-1; j++){
            B(i,j) = 0;
        }
    }

    //Grensebetingelser
    for(int i =0; i<n; i++){
        A(0,i) = 0.0;
        A(n-1,i) = 0.0;
        A(i,0) = 0;

    }
    for(int i =0; i<n; i++){
         A(i,n-1) = 1.0;
    }
   // B.print();

    ofstream gfile;
    string gfilename = "RESULTS/two_dimension_evolution_of_time.dat";
    gfile.open(gfilename);
    gfile << setiosflags(ios::showpoint | ios::uppercase);
    /*
    for( int j = 0; j<n; j++){
        for(int i=0; i<n; i++){
         gfile << setw(20) << setprecision(8) << A(j,i);
         }
    }
    gfile<<endl;
    */
    //Starter iterative løser
    for(int k=0; k<MaxIterations; k++){
        //B.print();
        //cout << k << endl;

        for(int i = 1; i<n-1; i++){
            for(int j=1; j<n-1; j++){
                A(i,j) = B(i,j) + D*(B(i+1,j) + B(i,j+1) -4*B(i,j) + B(i-1,j) + B(i,j-1));

                }
            }
        for( int j = 0; j<n; j++){
            for(int i=0; i<n; i++){
             gfile << setw(20) << setprecision(8) << A(j,i);
             }
        }
        gfile <<endl;
        //cout << A << endl;

        double sum = 0.0;

        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                sum += (B(i,j)-A(i,j))*(B(i,j)-A(i,j));
                B(i,j) = A(i,j);
            }
        }
        if(sqrt(sum)<abstol){
            gfile.close();
            return k;
        }
    }
    gfile.close();
    return MaxIterations;

}

void diffusjon2dim(int n){
    mat A =zeros<mat>(n,n);
    mat q = zeros<mat>(n,n);
    double tol = 10e-8;
    double dx = 1.0/(n-1);
    double dt = 0.01*dx*dx;
    cout << dx << endl;
    cout << dt << endl;
    /*for(int i =0; i < n; i++){
        for(int j = 0; j<n; j++){
            q(i,j) = 0;
        }
    }*/

    //for(int i =0; i<tsteps; i++){
        //cout << A << endl;
        int itcount = JacobiSolver(n,dt,dx,A,tol);
        cout << itcount << endl;
    //}
}


double Q_sim2(double x,double rho, double cp, double delta_t){
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
    return ((0.05+0.5)*pow(10,-6)*delta_t/(rho*cp)); // Need to be scaled appropriatly
  }
  cout<<"Check"<<endl;
 return 0;
}
void simulation_with_enrichment(int nx,int ny, double dt, double dxy, mat &A, int tsteps){
    clock_t start, finish;


    //Simulation spesific constants ;

    double rho = pow(3.510,3);      // Densety              ;       Kg/m^3
    double cp = 1000;               // Heat capasity        ;       J/Kg/ degree Celcius
    double k = 2.5;                 // Thermal conductivity ;       W/m/ degree Celcius


    //Backward Euler ;

    mat B = A; //Gammel A
    double abstol = 10e-4;
    double D = dt/(dxy*dxy)*k/(rho*cp);
    double scaled_x;
    ofstream gfile;
    string gfilename = "RESULTS/simulation_with_enrichment.dat";
    gfile.open(gfilename);
    gfile << setiosflags(ios::showpoint | ios::uppercase);

    for( int i = 0; i<nx; i++){
        for(int j=0; j<ny; j++){
         gfile << setw(20) << setprecision(8) << A(i,j);
         }
    }
    gfile<<endl;
    start = clock();
    //Starter iterative løser
    for(int k=0; k<tsteps; k++){

        for(int i = 1; i<nx-1; i++){
            scaled_x =(double)i/ (double)nx;
            for(int j=1; j<ny-1; j++){

                A(i,j) = ( (B(i,j) + D*(B(i+1,j) + B(i,j+1) -4*B(i,j) + B(i-1,j) + B(i,j-1))))+Q_sim2(scaled_x,rho,cp,dt);
            }
        }


        //cout << A << endl;

        double sum = 0.0;

        for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                sum += (B(i,j)-A(i,j))*(B(i,j)-A(i,j));
                //B(i,j) = A(i,j);
            }
        }
        if(sqrt(sum)<abstol){
            //gfile.close();
            //return k;
            finish = clock();
            cout<<"Equilibrium reached at"<<k*dt/60/60/24/365/pow(10,9)<<" Giga years. (Simulation time ; "<< ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)/60<<" Minutes)" <<endl;
            k=tsteps;
        }
        B=A;
    }
    for( int i = 0; i<nx; i++){
        for(int j=0; j<ny; j++){
         gfile << setw(20) << setprecision(8) << A(i,j);
         }
    }

    gfile <<endl;
    //B=A;
    gfile.close();
    //return MaxIterations;

}


double Q_radioactive_decay(double x,double rho, double cp, double delta_t, double radio_active_contribution){
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


    return ((0.05+(radio_active_contribution))*pow(10,-6)*delta_t/(rho*cp)); // Need to be scaled appropriatly
  }
  cout<<"Check"<<endl;
 return 0;
}
void simulation_implemented_enrichment_decay(int nx,int ny, double dt, double dxy, mat &A, int tsteps,int start_step){
    //Simulation spesific constants ;

    double rho = pow(3.510,3);      // Densety              ;       Kg/m^3
    double cp = 1000;               // Heat capasity        ;       J/Kg/ degree Celcius
    double k = 2.5;                 // Thermal conductivity ;       W/m/ degree Celcius

    //Simulation spesific operators
    double time;   // The spesefic simulated time.
    double U;      // Radio contribution
    double Th;     //
    double K;
    double scaled_x;

    //Backward Euler ;

    mat B = A; //Gammel A
    double D = dt/(dxy*dxy)*k/(rho*cp);


    //Starter iterative løser
    for(int k=start_step; k<tsteps; k++){

        // Precalculate radioactive contribution in the manlte.
        time = k*dt/60/60/24/365/pow(10,9);   // Time in Giga years
        U  = 0.2 *pow(0.5 , time/4.47);       //  Half life  ;  4.47  Giga years;
        Th = 0.2 *pow(0.5 , time/14.0);       //  Half life  ; 14.0   Giga years;
        K  = 0.1 *pow(0.5 , time/1.25);       //  Half life  ;  1.25  Giga years;

        for(int i = 1; i<nx-1; i++){
            scaled_x =(double)i/ (double)nx;
            for(int j=1; j<ny-1; j++){

                A(i,j) = ( (B(i,j) + D*(B(i+1,j) + B(i,j+1) -4*B(i,j) + B(i-1,j) + B(i,j-1))))+Q_radioactive_decay(scaled_x,rho,cp,dt,(U+Th+K));
            }
        }


        B=A;
    }


    //return MaxIterations;

}

#endif // PDE_TWO_DIMENSIONAL_H

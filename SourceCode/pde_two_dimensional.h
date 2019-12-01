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


void diffusjon2dim(int n,int tsteps){
    mat A =zeros<mat>(n,n);
    mat q = zeros<mat>(n,n);
    double tol = 10e-8;
    double dx = 1.0/(n-1);
    double dt = 0.01*dx*dx;

    /*for(int i =0; i < n; i++){
        for(int j = 0; j<n; j++){
            q(i,j) = 0;
        }
    }*/

    //for(int i =0; i<tsteps; i++){
        //cout << A << endl;
        int itcount = JacobiSolver(n,dt,dx,A,tol);
        //cout << A << endl;
    //}
}

int JacobiSolver(int n, double dt, double dx,mat &A, double abstol){
    int MaxIterations = 50;
    mat B = zeros<mat>(n,n); //Gammel A

    double D = dt/(dx*dx);

    //Initialbetingelse
    for(int i=1; i <n-1; i++){
        for(int j=1; j<n-1; j++){
            B(i,j) = 1;
        }
    }

    //Grensebetingelser
    for(int i =0; i<n; i++){
        A(0,i) = 0.0;
        A(n-1,i) = 0.0;
        A(i,0) = 0;
        A(i,n-1) = 1.0;
    }
   // B.print();

    ofstream gfile;
    string gfilename = "results/crank_nicolson_dt001.dat";
    gfile.open(ffilename);
    gfile << setiosflags(ios::showpoint | ios::uppercase);
    //Starter iterative løser

    for(int k=0; k<MaxIterations; k++){
        //B.print();
        //cout << "HALLO" << endl;

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
         cout << A << endl;

        double sum = 0.0;

        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                sum += (B(i,j)-A(i,j))*(B(i,j)-A(i,j));
                B(i,j) = A(i,j);
            }
        }
        if(sqrt(sum)<abstol){
            return k;
        }
    }
    return MaxIterations;

}
#endif // PDE_TWO_DIMENSIONAL_H

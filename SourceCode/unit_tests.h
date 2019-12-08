#ifndef UNIT_TESTS_H
#define UNIT_TESTS_H
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <pde_single_dimensional.h>
#include <pde_two_dimensional.h>

using namespace std;
using namespace arma;


void test_forward_euler(){
    int n = 10;
    double delta_x =  0.1;
    double delta_t = 0.01;
    int tsteps = (int)(1/delta_t);
    double L = 1;

    //Analytical solution
    mat expected(tsteps,n+1);
    double x;
    double t;
    for(int j = 0; j<tsteps;j++){
        for(int i = 0; i<n;i++){
          t = delta_t*j;
          x = delta_x *i;
          expected(j,i) = x/L + (pow(-1,n) *  2/(n*M_PI)*sin((n*M_PI*x)/L)*exp(( (-n)*(-n)*M_PI*M_PI*t )/(L*L) ));
        }
     }

    mat approximation= forward_euler(delta_x,delta_t);
    //cout<<approximation<<endl;
    for(int i = 0; i<n;i++){
    cout<<expected(tsteps-1,i) - approximation(tsteps-1,i) <<endl;
    }
}


void test_jacobi_solver(){
    int n = 10;
    int tsteps = 10000000;

    mat A =zeros<mat>(n,n);
    mat q = zeros<mat>(n,n);
    double tol = 10e-8;
    double dx = 1.0/(n-1);
    double dt = 0.01*dx*dx;

    int itcount = JacobiSolver(n,dt,dx,A,tol);

    // Test the Soultution
    double sum = 0;
    double expected;
    for(int i = 0; i<n;i++){
        for(int j = 0; j<n; j++){
            expected = sin(dx*i*M_PI)*sin(dx*j*M_PI)*exp(-M_PI*M_PI*(i*i*j*j)*itcount*dt);
            cout<<expected<<endl;
            sum += fabs(A(i,j) - expected);
            cout<<fabs(A(i,j)) <<endl;
        }
    }
cout <<sum<<endl;
cout<<A<<endl;

}


#endif // UNIT_TESTS_H

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




vec analytical_solution(double t){
    int n = 11;
    vec x = linspace<vec>(0,1,n);

    vec res = zeros<vec>(n);

    for(int k = 1; k<n; k++){
        for(int i = 0;i<n; i++){
            res(i) = res(i)+ 2*pow(-1,k)/(k*M_PI)*exp(-k*k*M_PI*M_PI*t)*sin(M_PI*k*x(i));
        }
    }
    for(int i = 0;i<n; i++){
        res(i) = res(i)+ x(i);
    }
   return(res);
}


void test_forward_euler(){
    //
    double tol=10e-3;

    //Analytical solution
    vec expected;

    // Setup numerical test solution.
    int n = 10;
    double delta_x =  0.1;
    double delta_t = 0.001;
    double lenght = 1;
    double time = 1;

    // Initial condition.
    vec u = zeros<vec>(n+1);
    u(n) = 1;

    mat approximation = forward_euler(delta_x,delta_t, lenght, time, u);

    //test at time 0.1
    double t = 0.1;
    expected = analytical_solution(t);
    double sum = 0.0;
    for(int i=0; i<=n; i++){
      sum += (expected(i)-approximation((t/delta_t),i))*(expected(i)-approximation((t/delta_t),i));
    }

   // Check if  good enough.
    // Test that our results seems to be okay.
    try{
        if (sqrt(sum)>tol){
            throw "Warning: The Forward Euler function is not as accurate as expected, there may be a serious problem that needs to be checked up upon.";
        }
    }
    catch (const char* msg){
        cerr << msg <<endl;
    }

   //test at time 1
   tol = 10e-6;
   t = 1;
   expected = analytical_solution(t);
   sum = 0.0;
   for(int i=0; i<=n; i++){
     sum += (expected(i)-approximation((t/delta_t),i))*(expected(i)-approximation((t/delta_t),i));
   }
//   cout<<sqrt(sum)<<endl;

   // Check if  good enough.

   try{
       if (sqrt(sum)>tol){
           throw "Warning: The Forward Euler function is not as accurate as expected, there may be a serious problem that needs to be checked up upon.";
       }
   }
   catch (const char* msg){
       cerr << msg <<endl;
   }

}


void test_backward_euler(){
    //
    double tol=10e-3;

    //Analytical solution
    vec expected;

    // Setup numerical test solution.
    int n = 10;
    double delta_x =  0.1;
    double delta_t = 0.001;
    double lenght = 1;
    double time = 1;

    // Initial condition.
    vec u = zeros<vec>(n+1);
    u(n) = 1;

    mat approximation = backward_euler(delta_x,delta_t, lenght, time, u);

    //test at time 0.1
    double t = 0.1;
    expected = analytical_solution(t);
    double sum = 0.0;
    for(int i=0; i<=n; i++){
      sum += (expected(i)-approximation((t/delta_t),i))*(expected(i)-approximation((t/delta_t),i));
    }
   //cout<<sqrt(sum)<<endl;
   // Check if  good enough.
    // Test that our results seems to be okay.
    try{
        if (sqrt(sum)>tol){
            throw "Warning: The Backward Euler function is not as accurate as expected, there may be a serious problem that needs to be checked up upon.";
        }
    }
    catch (const char* msg){
        cerr << msg <<endl;
    }

   //test at time 1
   tol = 10e-6;
   t = 1;
   expected = analytical_solution(t);
   sum = 0.0;
   for(int i=0; i<=n; i++){
     sum += (expected(i)-approximation((t/delta_t),i))*(expected(i)-approximation((t/delta_t),i));
   }
   //cout<<sqrt(sum)<<endl;

   // Check if  good enough.

   try{
       if (sqrt(sum)>tol){
           throw "Warning: The Backward Euler function is not as accurate as expected, there may be a serious problem that needs to be checked up upon.";
       }
   }
   catch (const char* msg){
       cerr << msg <<endl;
   }

}


void test_crank_nicolson(){
    //
    double tol=10e-3;

    //Analytical solution
    vec expected;

    // Setup numerical test solution.
    int n = 10;
    double delta_x =  0.1;
    double delta_t = 0.001;
    double lenght = 1;
    double time = 1;

    // Initial condition.
    vec u = zeros<vec>(n+1);
    u(n) = 1;

    mat approximation = crank_nicolson(delta_x,delta_t, lenght, time, u);

    //test at time 0.1
    double t = 0.1;
    expected = analytical_solution(t);
    double sum = 0.0;
    for(int i=0; i<=n; i++){
      sum += (expected(i)-approximation((t/delta_t),i))*(expected(i)-approximation((t/delta_t),i));
    }
   cout<<sqrt(sum)<<endl;
   // Check if  good enough.
    // Test that our results seems to be okay.
    try{
        if (sqrt(sum)>tol){
            throw "Warning: The Crank Nicolson function is not as accurate as expected, there may be a serious problem that needs to be checked up upon.";
        }
    }
    catch (const char* msg){
        cerr << msg <<endl;
    }

   //test at time 1
   tol = 10e-6;
   t = 1;
   expected = analytical_solution(t);
   sum = 0.0;
   for(int i=0; i<=n; i++){
     sum += (expected(i)-approximation((t/delta_t),i))*(expected(i)-approximation((t/delta_t),i));
   }
   cout<<sqrt(sum)<<endl;

   // Check if  good enough.

   try{
       if (sqrt(sum)>tol){
           throw "Warning: The Crank Nicolson function is not as accurate as expected, there may be a serious problem that needs to be checked up upon.";
       }
   }
   catch (const char* msg){
       cerr << msg <<endl;
   }

}




void test_jacobi_solver(){

    /*
    Test our implemention of the Jacobi solver function , with a known analytical solution.
    */

    //Setup inital conditions.
    int n = 10;
    double tol = 10e-8;
    double dx = 1.0/(n-1);
    double dt = 0.01*dx*dx;
    mat A =zeros<mat>(n,n);
    for(int i = 1; i<n-1;i++){
        for(int j = 1; j<n-1; j++){
            A(i,j) = 1;
        }
    }

    // Calculation on our test.
    int itcount = JacobiSolver(n,n,dt,dx,A,tol);

    // Find analytical solution and compare it to our caculations.
    mat expected(n,n);
    double sum = 0;
    for(int i = 0; i<n;i++){
        for(int j = 0; j<n; j++){
            expected(i,j) = sin(dx*i*M_PI)*sin(dx*j*M_PI)*exp(-M_PI*M_PI*(i*i+j*j)*itcount*dt);
            sum += fabs(fabs(A(i,j)) - fabs(expected(i,j)));
        }
    }

    tol = 10e-4; // Set tolerance of the total error between our analytical solution and our calculation.
    // Test that our results seems to be okay.
    try{
        if (sum > tol){
            throw "Warning: The Jacobi solver function is not as accurate as expected, there may be a serious problem that needs to be checked up upon.";
        }
    }
    catch (const char* msg){
        cerr << msg <<endl;
    }

}

void test_package_installed(){

}

#endif // UNIT_TESTS_H

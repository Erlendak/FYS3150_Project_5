#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <pde_single_dimensional.h>
#include <pde_two_dimensional.h>
#include <unit_tests.h>

using namespace std;


int main(){

    //test_jacobi_solver();
    //test_forward_euler();
/*
    double delta_x1 = 0.1;
    double delta_x2 = 0.01;
    double delta_t = 0.000005;

    //Forward Euler

    //Delta X ; 0.1
    ofstream afile;
    string afilename = "RESULTS/forward_euler_dt01.dat";
    afile.open(afilename);
    afile << setiosflags(ios::showpoint | ios::uppercase);
    afile << setw(20) << setprecision(8) << forward_euler(delta_x1,delta_t) <<endl;
    afile.close();

    //Delta X ; 0.01
    ofstream bfile;
    string bfilename = "RESULTS/forward_euler_dt001.dat";
    bfile.open(bfilename);
    bfile << setiosflags(ios::showpoint | ios::uppercase);
    bfile << setw(20) << setprecision(8) << forward_euler(delta_x2,delta_t) <<endl;
    bfile.close();



    //Backward Euler ;

    // Delta X ; 0.1
    int tsteps = (int)(1/delta_t)-1;
    double alpha =  delta_t/(delta_x1*delta_x1);
    ofstream cfile;
    string cfilename = "RESULTS/backward_euler_dt01.dat";
    cfile.open(cfilename);
    cfile << setiosflags(ios::showpoint | ios::uppercase);
    cfile << setw(20) << setprecision(8) << backward_euler((1/delta_x1)-1, tsteps, delta_x1, alpha) <<endl;
    cfile.close();


    //Delta X ; 0.01
    alpha =  delta_t/(delta_x2*delta_x2);
    ofstream dfile;
    string dfilename = "RESULTS/backward_euler_dt001.dat";
    dfile.open(dfilename);
    dfile << setiosflags(ios::showpoint | ios::uppercase);
    dfile << setw(20) << setprecision(8) << backward_euler((1/delta_x2)-1, tsteps, delta_x2, alpha) <<endl;
    dfile.close();



    // Crank Nicolson ;


    //Delta X ; 0.1
    tsteps = (int)(1/delta_t)-1;
    alpha =  delta_t/(delta_x1*delta_x1);
    cout<<alpha<<endl;
    ofstream efile;
    string efilename = "RESULTS/crank_nicolson_dt01.dat";
    efile.open(efilename);
    efile << setiosflags(ios::showpoint | ios::uppercase);
    efile << setw(20) << setprecision(8) <<  crank_nicolson((1/delta_x1)-1, tsteps, delta_x1, alpha) <<endl;
    efile.close();

    //Delta X ; 0.01
    alpha =  delta_t/(delta_x2*delta_x2);
    cout<<alpha<<endl;
    ofstream ffile;
    string ffilename = "RESULTS/crank_nicolson_dt001.dat";
    ffile.open(ffilename);
    ffile << setiosflags(ios::showpoint | ios::uppercase);
    ffile << setw(20) << setprecision(8) <<  crank_nicolson((1/delta_x2)-1, tsteps, delta_x2, alpha)<<endl;
    ffile.close();


    diffusjon2dim(10);

*/
    int n= 120;
    int tsteps = 100000000;
    double delta_x = 1000;
    double delta_t = 1000000000;
    double alpha = delta_t/(delta_x*delta_x);
    ofstream hfile;
    string hfilename = "RESULTS/simulation_no_enrichment.dat";
    hfile.open(hfilename);
    hfile << setiosflags(ios::showpoint | ios::uppercase);
    vec no_enrichment = simulation_before_radioactive_enrichment(n, tsteps, delta_x,delta_t);
    hfile << setw(20) << setprecision(8) << no_enrichment <<endl;
    hfile.close();

    double time = tsteps * delta_t/60/60/24/365;
    cout<< time/pow(10,9) <<endl;//in Giga years

    //Grensebetingelser
   int nx = 120;
   int ny = 150;
   double dxy=1000;
   double dt = delta_t;
   mat A = zeros<mat>(nx,ny);

   // Inital condtions
   for(int j =0; j<ny; j++){
       for(int i =0; i<nx; i++){
           A(i,j)    = no_enrichment(i); //Depth condtions use vector from task above
           //A(i,ny-1) = 1; //Depth condtions use vector from task above
        }
    }
   for(int j =0; j<ny; j++){
       A(0,j)    = 8; //Depth condtions use vector from task above
       A(nx-1,j) = 1300; //Depth condtions use vector from task above
    }
   tsteps = 1000000;


    simulation_with_enrichment(nx, ny, dt, dxy,A, tsteps);

    cout<<"Starting main simulation..." <<endl;

    // Main simulation


    ofstream ifile;
    string ifilename = "RESULTS/simulation_implemented_enrichment_decay.dat";
    ifile.open(ifilename);
    ifile << setiosflags(ios::showpoint | ios::uppercase);

    int n_results = 10;
    tsteps = pow(10,9)*365*24*60*60/dt/n_results; // Gives us amount of time steps we need to use per result we want to note of a total of a Giga year.
    int start_step;
    clock_t start, finish;
    for (int t = 0; t <n_results; t++){
        start = clock();
        start_step = t*tsteps;
        simulation_implemented_enrichment_decay(nx, ny, dt, dxy,A, (t+1)*tsteps,start_step);
        for( int i = 0; i<nx; i++){
            for(int j=0; j<ny; j++){
                ifile << setw(20) << setprecision(8) << A(i,j);
            }
        }
        ifile<<endl;
        finish = clock();
        cout<<"Simulated for ; "<<(t+1)*tsteps*dt/60/60/24/365/pow(10,9)<<" Giga years (Simulation time ; "<< ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)/60<<" Minutes)"  <<endl;
    }

    ifile.close();
    return 0;

};

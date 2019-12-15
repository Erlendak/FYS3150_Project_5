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

    test_jacobi_solver();
    test_forward_euler();
    test_backward_euler();
    test_crank_nicolson();

    //Initialize variabels used multiple times ;
    clock_t start, finish; // initialize time operators.

    // Setup examples ;
    double delta_x1 = 0.1;
    double delta_x2 = 0.01;
    double delta_t = 0.000005;
    double example_lenght = 1;
    double example_time = 1;
    int example_n1 =  (int)(example_lenght/delta_x1);
    int example_n2 = (int)(example_lenght/delta_x2);

    // Setup inital conditions for example.
    vec u1 = zeros<vec>(example_n1+1);     // For example with delta x = 0.1.
    vec u2 = zeros<vec>(example_n2+1);     // For example with delta x = 0.01.
    u1(example_n1) = u2(example_n2) = 1;   // Boundary condition.




    //Forward Euler


    //Delta X ; 0.1

    ofstream afile;
    string afilename = "RESULTS/forward_euler_dt01.dat";
    afile.open(afilename);
    afile << setiosflags(ios::showpoint | ios::uppercase);
    start = clock();
    afile << setw(20) << setprecision(8) << forward_euler(delta_x1,delta_t, example_lenght, example_time, u1) <<endl;
    finish = clock();
    afile.close();

    cout<<"Example simulation using Forward Euler with a time step of ; "<< delta_t << " seconds and a delta x ; "<< delta_x1<<
    " Meters.\n(Simulation time ; " << ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)<<" Seconds)\n" <<endl;


    //Delta X ; 0.01

    ofstream bfile;
    string bfilename = "RESULTS/forward_euler_dt001.dat";
    bfile.open(bfilename);
    bfile << setiosflags(ios::showpoint | ios::uppercase);
    start = clock();
    bfile << setw(20) << setprecision(8) << forward_euler(delta_x2,delta_t, example_lenght, example_time, u2) <<endl;
    finish = clock();
    bfile.close();

    cout<<"Example simulation using Forward Euler with a time step of ; "<< delta_t << " seconds and a delta x ; "<< delta_x2<<
    " Meters.\n(Simulation time ; " << ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)<<" Seconds)\n" <<endl;




    //Backward Euler ;


    // Delta X ; 0.1

    ofstream cfile;
    string cfilename = "RESULTS/backward_euler_dt01.dat";
    cfile.open(cfilename);
    cfile << setiosflags(ios::showpoint | ios::uppercase);
    start = clock();
    cfile << setw(20) << setprecision(8) << backward_euler(delta_x1,delta_t, example_lenght, example_time, u1) <<endl;
    finish = clock();
    cfile.close();

    cout<<"Example simulation using Backward Euler with a time step of ; "<< delta_t << " seconds and a delta x ; "<< delta_x1<<
    " Meters.\n(Simulation time ; " << ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)<<" Seconds)\n" <<endl;

    //Delta X ; 0.01

    ofstream dfile;
    string dfilename = "RESULTS/backward_euler_dt001.dat";
    dfile.open(dfilename);
    dfile << setiosflags(ios::showpoint | ios::uppercase);
    start = clock();
    dfile << setw(20) << setprecision(8) << backward_euler(delta_x2,delta_t, example_lenght, example_time, u2) <<endl;
    finish = clock();
    dfile.close();

    cout<<"Example simulation using Backward Euler with a time step of ; "<< delta_t << " seconds and a delta x ; "<< delta_x2<<
    " Meters.\n(Simulation time ; " << ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)<<" Seconds)\n" <<endl;




    // Crank Nicolson ;


    //Delta X ; 0.1

    ofstream efile;
    string efilename = "RESULTS/crank_nicolson_dt01.dat";
    efile.open(efilename);
    efile << setiosflags(ios::showpoint | ios::uppercase);
    start = clock();
    efile << setw(20) << setprecision(8) <<  crank_nicolson(delta_x1,delta_t, example_lenght, example_time, u1) <<endl;
    finish = clock();
    efile.close();

    cout<<"Example simulation using Crank Nicolson with a time step of ; "<< delta_t << " seconds and a delta x ; "<< delta_x1<<
    " Meters.\n(Simulation time ; " << ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)<<" Seconds)\n" <<endl;

    //Delta X ; 0.01

    ofstream ffile;
    string ffilename = "RESULTS/crank_nicolson_dt001.dat";
    ffile.open(ffilename);
    ffile << setiosflags(ios::showpoint | ios::uppercase);
    start = clock();
    ffile << setw(20) << setprecision(8) <<  crank_nicolson(delta_x2,delta_t, example_lenght, example_time, u2)<<endl;
    finish = clock();
    ffile.close();

    cout<<"Example simulation using Crank Nicolson with a time step of ; "<< delta_t << " seconds and a delta x ; "<< delta_x2<<
    " Meters.\n(Simulation time ; " << ( ( (double)finish - (double)start ) /CLOCKS_PER_SEC)<<" Seconds)\n" <<endl;

/*
    diffusjon_example_2dim(10);


    int n= 120;
    int tsteps = 100000000;
    double delta_x = 1000;
    double delta_t = 1000000000;
    double alpha = delta_t/(delta_x*delta_x);
    ofstream hfile;
    string hfilename = "RESULTS/simulation_no_enrichment.dat";
    hfile.open(hfilename);
    hfile << setiosflags(ios::showpoint | ios::uppercase);
    vec no_enrichment = simulation_no_radioactive_enrichment(n, tsteps, delta_x,delta_t);
    hfile << setw(20) << setprecision(8) << no_enrichment <<endl;
    hfile.close();

    //Grensebetingelser
   int nx = 120;
   int ny = 150;
   double dxy=1000;
   double dt = delta_t;
   mat A = zeros<mat>(nx,ny);

   // Inital conditions
   for(int j =0; j<ny; j++){
       for(int i =0; i<nx; i++){
           A(i,j)    = no_enrichment(i); //Depth conditions use vector from task above
           //A(i,ny-1) = 1; //Depth conditions use vector from task above
        }
    }
   for(int j =0; j<ny; j++){
       A(0,j)    = 8; //Depth conditions use vector from task above
       A(nx-1,j) = 1300; //Depth conditions use vector from task above
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
*/
    return 0;

};

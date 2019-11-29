#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>
#include <pde_single_dimensional.h>

using namespace std;


int main(){
    double delta_x1 = 0.1;
    double delta_x2 = 0.01;
    double delta_t = 0.001;
    //forward_euler(delta_x1);
    //forward_euler(delta_x2);

    //forward_euler(delta_x1,0.001);
    int tsteps = (int)(1/delta_t)-1;
    double alpha =  delta_t/(delta_x1*delta_x1);
    backward_euler((1/delta_x1)-1, tsteps, delta_x1, alpha);
    //crank_nicolson(1/delta_x1, tsteps, delta_x1, alpha);

    return 0;

};

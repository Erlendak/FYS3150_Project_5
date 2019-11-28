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

    //forward_euler(delta_x1);
    //forward_euler(delta_x2);
    barbarian_forward_euler(delta_x1,0.001);


    return 0;

};

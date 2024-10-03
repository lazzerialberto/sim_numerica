#include "FunzioneBase.h"
#include "utils.h"
#include "random.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

using namespace std;


int main(int argc, char* argv[]){

    if(argc!=4){
        cout << "ERROR: \nStarting Temperature: \nDelta mu: \nDelta sigma:" << endl;
        exit(-1);
    }

    int N=100;
    int M=500000;
    int SA_steps=1000;

    double mu=1.;
    double sigma=1.;
    double init_temp=atof(argv[1]);
    double delta_mu=atof(argv[2]);
    double delta_sigma=atof(argv[3]);
    Random rnd;
    rnd.Initialize();
    Boltzmann boltz;
    int count=0;

    SAnnealing simulation(N,M,mu,sigma,init_temp,delta_mu,delta_sigma,rnd);

    while(simulation.get_T()>=0.001){
        simulation.Step(count,SA_steps,boltz);
        count++;
    }

    simulation.finalize(count);

    rnd.SaveSeed();
    return 0;
}
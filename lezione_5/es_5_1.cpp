#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <memory>

#include "random.h"
#include "utils.h"
#include "FunzioneBase.h"

using namespace std;

int main(int argc, char * argv[]){

    if(argc!=8){
        std::cout << "Errors:\ntransition probability needed, write gauss or unif\nset transition pass\nset ground (g) or excited (e) state\nset starting position\nequilibration nedeed? (yes or no)" << endl;
        exit(-1);
    }

    int N=100;
    int M=100000;
    bool gauss;
    groundwf ground;
    excitedwf excited;

    double delta=atof(argv[2]);
    posizione r(atof(argv[4]),atof(argv[5]),atof(argv[6]));

    FunzioneBase &prob= (std::string(argv[3])=="g") ? static_cast<FunzioneBase&>(ground) : static_cast<FunzioneBase&>(excited);

    if(std::string(argv[1])=="gauss"){
        gauss=true;
    }
    else if(std::string(argv[1])=="unif"){
        gauss=false;
    }

    Random rnd;
    rnd.Initialize();

    Metropolis metro(delta,gauss,rnd);

    //setting delta in order to have ~50% acceptance
    metro.SetStepLenght(r,prob);

    posizione r_n(atof(argv[4]),atof(argv[5]),atof(argv[6]));

    blockingaverage block(argv,M);

    //equilibration
    if(std::string(argv[7])=="yes"){
        for(int j=0; j<40; j++){
            block.measure(M,metro,prob,r_n);
        }
    }

    for(int i=0;i<N; i++){
        block.averages(i,M,metro,prob,r_n);
    }

    rnd.SaveSeed();

    return 0;
}
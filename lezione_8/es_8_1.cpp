#include <iostream>
#include <iomanip>
#include <fstream>

#include "utils.h"
#include "random.h"
#include "FunzioneBase.h"

using namespace std;

void Progress_Bar(int progress, int total, int bar_width=50) {

    float percentage = static_cast<float>(progress) / total;
    int pos = static_cast<int>(bar_width * percentage);

    std::string bar;
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) bar += "=";
        else if(i==pos) bar += ">";
        else bar += " ";
    }

    std::cout<< "|" << bar << "| " << setprecision(3) << int(percentage * 100.0)<< "%\r";
    std::fflush(stdout);
}

int main(int argc, char * argv[]){

    if(argc!=6){
        std::cout << "Errors:\ntransition probability needed, write gauss or unif\nset transition pass\nset starting position\nset mu\nset sigma" << endl;
        exit(-1);
    }

    int N=100;
    int M=50000;
    bool gauss;
    double delta=atof(argv[2]);
    double x_start=atof(argv[3]);
    double x=atof(argv[3]);
    double mu=atof(argv[4]);
    double sigma=atof(argv[5]);

    gausswf prob;

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
    metro.SetStepLenght(x,prob,mu,sigma);

    blockingaverage block(argv,M);

    for(int i=0;i<N; i++){
        block.averages(true,i,M,metro,prob,x,mu,sigma);
        Progress_Bar(i,N);
    }

    cout << endl;
    block.measure_savepos(M,metro,prob,x_start,mu,sigma);

    rnd.SaveSeed();

    return 0;
}
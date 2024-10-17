#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#include "utils.h"
#include "random.h"
#include "posizione.h"

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


int main(){

    Random rnd;
    rnd.Initialize();

    TSP tsp;

    Population pop(tsp.get_individuals(),tsp.get_cities());
    pop.initial_pop(rnd,tsp.get_circle());
    pop.Order_pop();
    pop.print_best_path(1);
    pop.print_Losses(1);

    for(int i=1; i<tsp.get_generations(); i++){
        pop.selection_and_crossover(rnd,tsp.get_p_sel(), tsp.get_p_cross());
        pop.Mutation(rnd,tsp.get_p_pp(), tsp.get_p_s(), tsp.get_p_p(),tsp.get_p_i());
        pop.Order_pop();
        pop.print_Losses(i+1);
        pop.print_best_path(i+1);
        Progress_Bar(i,tsp.get_generations());
    }
    std::cout << endl;

    rnd.SaveSeed();
    return 0;
}
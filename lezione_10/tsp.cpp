#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

#include <mpi.h>

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


int main(int argc,char*argv[]){

    Random rnd;
    rnd.Initialize(0);
    TSP tsp;
    element e(tsp.get_cities());
    e.generate_initial_config(rnd,tsp.get_circle());
    int n_migr=1; //number of element to exchange in the migration

    MPI_Init(&argc,&argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Ottieni il rank del processo corrente
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    double t_start= MPI_Wtime();
    
    rnd.Initialize(world_rank);

    tsp.print_output(world_rank,world_size);

    Population pop(tsp.get_individuals(),tsp.get_cities());
    pop.initial_pop_mpi(rnd,tsp.get_circle(),e);
    pop.Order_pop();
    if(world_rank==0){
        pop.print_best_path(1);
        pop.print_Losses(1);
    }

    for(int i=1; i<tsp.get_generations(); i++){
        pop.selection_and_crossover(rnd,tsp.get_p_sel(), tsp.get_p_cross());
        pop.Mutation(rnd,tsp.get_p_pp(), tsp.get_p_s(), tsp.get_p_p(),tsp.get_p_i());
        pop.Order_pop();
        if(world_rank==0){
            pop.print_Losses(i+1);
            pop.print_best_path(i+1);
        }

        if(i%50==0){ //every 50 generations send the first n_migr configurations at the following process, and recive the first n_migr from the previous process

            int prev = (world_rank - 1 + world_size) % world_size;
            int next = (world_rank + 1) % world_size;

            for(int j=0; j<n_migr; j++){

                std::vector<int> message_to_send=pop.get_element(j).get_label();
                std::vector<int> message_received(tsp.get_cities());

                MPI_Request send_request, recv_request;
                MPI_Status status;

                MPI_Irecv(message_received.data(), message_received.size(), MPI_INT, prev, j, MPI_COMM_WORLD, &recv_request);
                MPI_Isend(message_to_send.data(), message_to_send.size(), MPI_INT, next, j, MPI_COMM_WORLD, &send_request);

                MPI_Wait(&recv_request, &status);
                MPI_Wait(&send_request, &status);

                element recieved(pop.get_element(j));
                recieved.set_label(message_received);

                pop.Set_element(j,recieved);
            }
        }

        //Progress_Bar(i,tsp.get_generations());
    }
    //std::cout << endl;

    double t_end = MPI_Wtime();
    double dt = t_end-t_start;

    tsp.finalize(world_rank,dt);

    MPI_Finalize();

    rnd.SaveSeed();
    return 0;
}
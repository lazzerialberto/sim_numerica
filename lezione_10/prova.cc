#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

    int my_values[3];
    for(int i=0;i<3;i++){
        my_values[i]=i+1;
    }
    // Inizializza l'ambiente MPI
    MPI_Init(&argc, &argv);

    // Ottieni il numero totale di processi
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Ottieni il rank del processo corrente
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Ottieni il nome del processore
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Stampa un messaggio da ogni processo
    std::cout << "Process " << world_rank << " out of " << world_size
              << " on " << processor_name << std::endl;


    cout<< "Prima: "<< my_values[0]<< " "<< my_values[1]<<" "<< my_values[2]<< " per il processo "<< world_rank<< endl;
    MPI_Bcast(my_values,3,MPI_INTEGER,0, MPI_COMM_WORLD);
    cout<< "Dopo: "<< my_values[0]<< " "<< my_values[1]<< " "<< my_values[2]<< " per il processo "<< world_rank<< endl;

    // Termina l'ambiente MPI
    MPI_Finalize();
    return 0;
}


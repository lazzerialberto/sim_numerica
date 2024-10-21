#include <mpi.h>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Ottiene il rank del processo
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Ottiene il numero di processi

    if (size != 4) {
        std::cerr << "Questo programma richiede esattamente 4 processi." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int prev = (rank - 1 + size) % size;  // Rank del processo precedente
    int next = (rank + 1) % size;         // Rank del processo successivo

    std::vector<int> message_to_send(5, rank);  // Un vettore con 5 elementi, tutti uguali al rank del processo
    std::vector<int> message_received(5);       // Un vettore per ricevere 5 interi

    MPI_Request send_request, recv_request;
    MPI_Status status;

    // Ogni processo riceve un vettore di interi dal rank precedente
    MPI_Irecv(message_received.data(), message_received.size(), MPI_INT, prev, 0, MPI_COMM_WORLD, &recv_request);

    // Ogni processo invia un vettore di interi al rank successivo
    MPI_Isend(message_to_send.data(), message_to_send.size(), MPI_INT, next, 0, MPI_COMM_WORLD, &send_request);

    // Attende la ricezione del vettore
    MPI_Wait(&recv_request, &status);
    std::cout << "Processo " << rank << " ha ricevuto il vettore da " << prev << ": ";
    for (int val : message_received) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    // Attende il completamento dell'invio
    MPI_Wait(&send_request, &status);

    MPI_Finalize();
    return 0;
}

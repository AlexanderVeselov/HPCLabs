/*

3. –еализовать рассылку значени€ n процессам с помощью двухточеч-
ных обменов. Ёффективность реализации сравнить с функцией MPI_Bcast().

*/

#include <mpi.h>
#include <iostream>

void Broadcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
    int comm_size;
    MPI_Comm_size(comm, &comm_size);
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    if (my_rank == root)
    {
        for (int i = 0; i < comm_size; ++i)
        {
            if (i == root) continue;
            MPI_Send(buffer, count, datatype, i, 0, comm);
        }
    }
    else
    {
        MPI_Status status;
        MPI_Recv(buffer, count, datatype, root, 0, comm, &status);
    }
}

void Test_MPI_Bcast()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (rank == 0)
    {
        std::cout << "MPI_Bcast" << std::endl;
    }

    double start_time = -1;
    if (rank == 0)
    {
        start_time = MPI_Wtime();
    }

    MPI_Bcast(&start_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank != 0)
    {
        std::cout << "Node " << rank << " received a message after: " << MPI_Wtime() - start_time << std::endl;
    }
}

void Test_Custom_Bcast()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    if (rank == 0)
    {
        std::cout << "Custom Bcast" << std::endl;
    }

    double start_time = -1;
    if (rank == 0)
    {
        start_time = MPI_Wtime();
    }

    Broadcast(&start_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank != 0)
    {
        std::cout << "Node " << rank << " received a message after: " << MPI_Wtime() - start_time << std::endl;
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    Test_MPI_Bcast();

    MPI_Barrier(MPI_COMM_WORLD);

    Test_Custom_Bcast();

    MPI_Finalize();

    return 0;
}

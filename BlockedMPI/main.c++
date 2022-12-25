#include "Blocked.h"
#include "Core.h"

#include <mpi.h>
#include <unistd.h>

#include <vector>

using namespace MPIcore;

int main(int argc, char *argv[])
{
    /* int key{1};
       while (key == 1)
       {
           sleep(5);
       } */

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &Core::MPI_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &Core::MPI_rank);

    size_t Bnum{Core::MPI_size};
    size_t Bsize{static_cast<size_t>(atoi(argv[1]))};

    Blocked blckd{Bsize, Bnum};

    if (Core::MPI_rank == 0)
        blckd.fill();

    blckd.findLU();

    if (Core::MPI_rank == (Core::MPI_size - 1))
    {
        blckd.check();

        blckd.findDet();

        blckd.printResult();
    }

    MPI_Finalize();

    return 0;
}

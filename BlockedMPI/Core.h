#ifndef CORE_H
#define CORE_H

#include <vector>
#include <iostream>
#include <chrono>

namespace MPIcore
{
    using vector = std::vector<double>;
    using vector2 = std::vector<vector>;
    using vector3 = std::vector<vector2>;
    using vector4 = std::vector<vector3>;

    using size_t = std::size_t;
    using std::cout, std::endl;

    using std::chrono::high_resolution_clock,
        std::chrono::duration_cast,
        std::chrono::duration,
        std::chrono::milliseconds;

    class Core
    {
    public:
        static int MPI_size;
        static int MPI_rank;

        // Init MPI
        static void init(int *argc, char ***argv);

        // Print matrix in terminal
        static void print(vector4 const &_Bmatrix, size_t _Bsize, size_t _size, std::string const _label);

        // Transpose matrix
        static void transpose(vector const &_matrix, vector &_forResult, size_t _size);

        //
        static void getRow(vector const &_matrix, vector &_forResult, size_t _size, size_t _index);
        static void getRowAll(vector const &_matrix, vector &_forResult, size_t _size, int _rank, int rowsPerProc);

        static void getCol(vector const &_matrix, vector &_forResult, size_t _size, size_t _index);
        static void getColAll(vector const &_matrix, vector &_forResult, size_t _size, int _rank, int rowsPerProc);

        static void fillRowRank(vector const &_matrix, vector &forResult, size_t _size, int rank);
    };
}
#endif // CORE_H
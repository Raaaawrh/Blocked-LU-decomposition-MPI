#include "Core.h"

#include <iomanip>

using namespace MPIcore;

int Core::MPI_rank{0};
int Core::MPI_size{0};

void Core::init(int *argc, char ***argv)
{
}

// Show Matrix in Terminal
void Core::print(vector4 const &_Bmatrix, size_t _Bsize, size_t _size, std::string const _label)
{
    cout << _label << endl;
    for (size_t row{0}; row < _size; row++)
    {
        for (size_t col{0}; col < _size; col++)
        {
            cout << _Bmatrix[row / _Bsize][col / _Bsize][row % _Bsize][col % _Bsize] << '\t';
        }
        cout << endl;
    }
    cout << endl;
}

// Transpose Matrix
void Core::transpose(vector const &_matrix, vector &_forResult, size_t _size)
{
    for (int row = 0; row < _size; row++)
        for (int col = 0; col < _size; col++)
            _forResult[row * _size + col] = _matrix[col * _size + row];
}

// Get Row from Matrix
void Core::getRow(vector const &_matrix, vector &_forResult, size_t _size, size_t _index)
{
    for (int i = 0; i < _size; i++)
        _forResult[i] = _matrix[_index * _size + i];
}

// Get All Rows for Process
void Core::getRowAll(vector const &_matrix, vector &_forResult, size_t _size, int _rank, int _rowsPerProc)
{
    for (int sub_row = 0; sub_row < _rowsPerProc; sub_row++)
        for (int col_index = 0; col_index < _size; col_index++)
            _forResult[sub_row * _size + col_index] = _matrix[(_rank * _rowsPerProc * _size) + (sub_row * _size) + col_index];
}

// Get Column from Matrix
void Core::getCol(vector const &_matrix, vector &_forResult, size_t _size, size_t _index)
{
    for (int i = 0; i < _size; i++)
    {
        _forResult[i] = _matrix[_index * _size + i];
    }
}

// Get All Columns for Process
void Core::getColAll(vector const &_matrix, vector &_forResult, size_t _size, int _rank, int _colsPerProc)
{
    int index;
    for (int sub_col = 0; sub_col < _colsPerProc; sub_col++)
    {
        for (int row_index = 0; row_index < _size; row_index++)
        {
            _forResult[sub_col * _size + row_index] = _matrix[(row_index * _size) + (_rank * _colsPerProc) + sub_col];
        }
    }
}

// Filling Rows for Processes
void Core::fillRowRank(vector const &_matrix, vector &forResult, size_t _size, int _rank)
{
    for (int i = 0; i < _size; i++)
    {
        forResult[i] = _matrix[_rank * _size + i];
    }
}
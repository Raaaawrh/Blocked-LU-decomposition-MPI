#include "Blocked.h"

#include <random>
#include <mpi.h>

Blocked::Blocked(size_t _Bsize, size_t _Bnum)
    : m_size{_Bsize * _Bnum},
      m_Bsize{_Bsize}, m_Bnum{_Bnum},
      m_Bmatrix{vector4(m_Bnum,
                        vector3(m_Bnum,
                                vector2(m_Bsize,
                                        vector(m_Bsize,
                                               0.0))))},
      m_BmatrixCheck{vector4(m_Bnum,
                             vector3(m_Bnum,
                                     vector2(m_Bsize,
                                             vector(m_Bsize,
                                                    0.0))))},
      m_inBuffer{vector(m_Bsize * m_Bsize, 0.0)},
      m_outBuffer{vector(m_Bsize * m_Bsize, 0.0)}
{
    if (Core::MPI_rank != (Core::MPI_size - 1))
    {
        for (vector3 &Brow : m_BmatrixCheck)
        {
            for (vector2 &Block : Brow)
            {
                for (vector &row : Block)
                    row.clear();
                Block.clear();
            }
            Brow.clear();
        }
        m_BmatrixCheck.clear();
    }
}

void Blocked::fill()
{
    std::uniform_real_distribution<double> distribDiag(0.9, 1.1);
    std::random_device rd;
    std::mt19937 mersenne(rd());

    for (size_t row{0}; row < m_size; row++)
    {
        for (size_t col{0}; col < row; col++)
            m_Bmatrix[row / m_Bsize][col / m_Bsize][row % m_Bsize][col % m_Bsize] = 1.0 / (row + 1);

        m_Bmatrix[row / m_Bsize][row / m_Bsize][row % m_Bsize][row % m_Bsize] = distribDiag(mersenne);

        for (size_t col{row + 1}; col < m_size; col++)
            m_Bmatrix[row / m_Bsize][col / m_Bsize][row % m_Bsize][col % m_Bsize] = 1.0 / (row + 1);
    }
}

void Blocked::printMatrix()
{
    Core::print(m_Bmatrix, m_Bsize, m_size, "Matrix");
}

void Blocked::luD(size_t _Brow)
{
    vector2 &block{m_Bmatrix[_Brow][_Brow]};

    for (size_t diag{0}; diag < m_Bsize; diag++)
    {
        for (size_t row{diag + 1}; row < m_Bsize; row++)
        {
            block[row][diag] /= block[diag][diag];

            for (size_t col{diag + 1}; col < m_Bsize; col++)
                block[row][col] -= block[row][diag] * block[diag][col];
        }
    }
}

void Blocked::luR(size_t _Brow, size_t _Bcol)
{
    vector2 &block{m_Bmatrix[_Brow][_Bcol]};
    for (size_t diag{0}; diag < m_Bsize; diag++)
    {
        for (size_t row{diag + 1}; row < m_Bsize; row++)
            for (size_t col{0}; col < m_Bsize; col++)
                block[row][col] -= m_Bmatrix[_Brow][_Brow][row][diag] * block[diag][col];
    }
}

void Blocked::luC(size_t _Brow, size_t _Bcol)
{
    vector2 &block{m_Bmatrix[_Brow][_Bcol]};
    for (size_t diag{0}; diag < m_Bsize; diag++)
    {
        for (size_t row{0}; row < m_Bsize; row++)
        {
            block[row][diag] /= m_Bmatrix[_Bcol][_Bcol][diag][diag];
            for (size_t col{diag + 1}; col < m_Bsize; col++)
                block[row][col] -= m_Bmatrix[_Bcol][_Bcol][diag][col] * block[row][diag];
        }
    }
}

void Blocked::luI(size_t _Bdiag, size_t _Brow, size_t _Bcol)
{
    vector2 &block{m_Bmatrix[_Brow][_Bcol]};

    for (size_t diag{0}; diag < m_Bsize; diag++)
    {
        for (size_t row{0}; row < m_Bsize; row++)
            for (size_t col{0}; col < m_Bsize; col++)
                block[row][col] -= m_Bmatrix[_Brow][_Bdiag][row][diag] * m_Bmatrix[_Bdiag][_Bcol][diag][col];
    }
}

void Blocked::findLU()
{
    MPI_Barrier(MPI_COMM_WORLD);

    // Timer
    if (Core::MPI_rank == (Core::MPI_size - 1))
    {
        m_timeStart = MPI_Wtime();
    }

    // Bcast initial matrix
    if (Core::MPI_rank == 0)
    {
        for (int row = 1; row < m_Bnum; row++)
            for (int col = 0; col < m_Bnum; col++)
                sendBlock(row, col, row);
    }
    else
    {
        for (int col = 0; col < m_Bnum; col++)
            recvBlock(Core::MPI_rank, col, 0);
    }

    // Start LU
    for (size_t Bdiag{0}; Bdiag < m_Bnum; Bdiag++)
    {
        // Diag Block
        if (Core::MPI_rank == static_cast<int>(Bdiag))
        {
            luD(Bdiag);

            // Send Diag Block to All
            for (int rank{Bdiag + 1}; rank < m_Bnum; rank++)
                sendBlock(Bdiag, Bdiag, rank);

            // Row decomposition
            for (size_t Bcol{Bdiag + 1}; Bcol < m_Bnum; Bcol++)
            {
                luR(Bdiag, Bcol);
                // Send Row Blocks to All
                for (int rank{Bdiag + 1}; rank < m_Bnum; rank++)
                    sendBlock(Bdiag, Bcol, rank);
            }
        }
        else if (Core::MPI_rank > Bdiag)
        {
            // Get Diagonal Block
            recvBlock(Bdiag, Bdiag, Bdiag);
            luC(Core::MPI_rank, Bdiag);

            for (size_t BsubCol{Bdiag + 1}; BsubCol < m_Bnum; BsubCol++)
            {
                recvBlock(Bdiag, BsubCol, Bdiag);
                luI(Bdiag, Core::MPI_rank, BsubCol);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (Core::MPI_rank == (Core::MPI_size - 1))
    {
        m_timeStop = MPI_Wtime();

        m_timeLU = m_timeStop - m_timeStart;
    }
}

void Blocked::check()
{
    for (int row{0}; row < m_size; row++)
    {
        for (int col{0}; col < m_size; col++)
        {
            int count{std::min(row - 1, col)};
            for (int i{0}; i <= count; i++)
                m_BmatrixCheck[row / m_Bsize][col / m_Bsize][row % m_Bsize][col % m_Bsize] +=
                    m_Bmatrix[row / m_Bsize][i / m_Bsize][row % m_Bsize][i % m_Bsize] *
                    m_Bmatrix[i / m_Bsize][col / m_Bsize][i % m_Bsize][col % m_Bsize];
        }

        for (int col{row}; col < m_size; col++)
            m_BmatrixCheck[row / m_Bsize][col / m_Bsize][row % m_Bsize][col % m_Bsize] +=
                m_Bmatrix[row / m_Bsize][col / m_Bsize][row % m_Bsize][col % m_Bsize];
    }

    normL2();
}

void Blocked::normL2()
{
    double norm{0};

    for (size_t row{0}; row < m_size; row++)
        for (size_t col{0}; col < m_size; col++)
            norm += (m_Bmatrix[row / m_Bsize][col / m_Bsize][row % m_Bsize][col % m_Bsize] -
                     m_BmatrixCheck[row / m_Bsize][col / m_Bsize][row % m_Bsize][col % m_Bsize]) *
                    (m_Bmatrix[row / m_Bsize][col / m_Bsize][row % m_Bsize][col % m_Bsize] -
                     m_BmatrixCheck[row / m_Bsize][col / m_Bsize][row % m_Bsize][col % m_Bsize]) /
                    m_size / m_size;

    m_residual = std::sqrt(norm);
}

void Blocked::findDet()
{
    auto timeStart = high_resolution_clock::now();

    m_determinant = 1.0;
    for (size_t row{0}; row < m_size; row++)
        m_determinant *= m_Bmatrix[row / m_Bsize][row / m_Bsize][row % m_Bsize][row % m_Bsize];

    auto timeStop = high_resolution_clock::now();

    m_timeDet = duration<double, std::milli>(
                    duration_cast<milliseconds>(timeStop - timeStart))
                    .count();
}

void Blocked::sendBlock(size_t _Brow, size_t _Bcol, int const _destRank)
{
    for (size_t row{0}; row < m_Bsize; row++)
        for (size_t col{0}; col < m_Bsize; col++)
            m_outBuffer[row * m_Bsize + col] = m_Bmatrix[_Brow][_Bcol][row][col];

    MPI_Send(m_outBuffer.data(), m_outBuffer.size(), MPI_DOUBLE, _destRank, _destRank + 10, MPI_COMM_WORLD);
}

void Blocked::recvBlock(size_t _Brow, size_t _Bcol, int const _srceRank)
{
    MPI_Recv(m_inBuffer.data(), m_inBuffer.size(), MPI_DOUBLE, _srceRank, Core::MPI_rank + 10, MPI_COMM_WORLD, &status);

    for (size_t row{0}; row < m_Bsize; row++)
        for (size_t col{0}; col < m_Bsize; col++)
            m_Bmatrix[_Brow][_Bcol][row][col] = m_inBuffer[row * m_Bsize + col];
}

void Blocked::printResult()
{
    cout << m_timeLU << '\t' << m_residual << '\t' << m_determinant << endl;
}

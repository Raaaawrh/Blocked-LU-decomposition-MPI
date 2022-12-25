#include "Blocked.h"

#include <random>

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
                                                    0.0))))} {}

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

void Blocked::print()
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
    auto timeStart{high_resolution_clock::now()};

    for (size_t Bdiag{0}; Bdiag < m_Bnum; Bdiag++)
    {
        luD(Bdiag);

        for (size_t Bcol{Bdiag + 1}; Bcol < m_Bnum; Bcol++)
            luR(Bdiag, Bcol);

        for (size_t Brow{Bdiag + 1}; Brow < m_Bnum; Brow++)
        {
            luC(Brow, Bdiag);

            for (size_t BsubCol{Bdiag + 1}; BsubCol < m_Bnum; BsubCol++)
                luI(Bdiag, Brow, BsubCol);
        }
    }

    auto timeStop{high_resolution_clock::now()};

    m_timeLU = duration<double, std::milli>(
                   duration_cast<milliseconds>(
                       timeStop - timeStart))
                   .count();
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
        duration_cast<milliseconds>(timeStop - timeStart)).count();
}

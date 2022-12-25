#if !defined(BLOCKED_H)
#define BLOCKED_H

#include "Core.h"

using namespace MPIcore;

class Blocked
{
public:
    Blocked(size_t _Bsize, size_t _Bnum);
    ~Blocked(){};

    void fill();
    // void Bfill();
    void print();

    void findLU();
    void check();
    
    void findDet();

private:
    size_t m_size;
    size_t m_Bsize;
    size_t m_Bnum;

    vector4 m_Bmatrix;
    vector4 m_BmatrixCheck;

    double m_residual;
    double m_determinant;
    double m_timeLU;
    double m_timeDet;

    void luD(size_t _Bdiag);
    void luR(size_t _Brow, size_t _Bcol);
    void luC(size_t _Brow, size_t _Bcol);
    void luI(size_t _Bdiag, size_t _Brow, size_t _Bcol);

    void normL2();
};

#endif // BLOCKED_H

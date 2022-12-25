#include "Blocked.h"
#include "Core.h"

#include <vector>

using namespace MPIcore;

int main(int argc, char *argv[])
{
    size_t Bsize{4};
    size_t Bnum{200};

    Blocked blckd{Bsize, Bnum};
    blckd.fill();

    blckd.findLU();

    //blckd.check();

    blckd.findDet();

    return 0;
}




#ifndef LGMC_H
#define LGMC_H

//#include "mpi.h"
#include "stdio.h"
#include "memory.h"
#include "lattice.h"


namespace LGMC_NS {

class LGMC {
 public:
    Memory *memory;      // memory allocation functions
    Lattice *lattice;    // lattice for SEI


    LGMC(int argc, char ** argv);
    ~LGMC();
    void print_usage();

};

}

#endif

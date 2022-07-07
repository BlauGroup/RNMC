
#ifndef LGMC_H
#define LGMC_H

//#include "mpi.h"
#include "stdio.h"
#include "memory.h"
#include "lattice.h"
#include "../core/dispatcher.h"
#include "../GMC/sql_types.h"
#include "../GMC/reaction_network.h"

namespace LGMC_NS {

class LGMC {
 public:
    Memory *memory;      // memory allocation functions
    Lattice *lattice;    // lattice for SEI
    
    Dispatcher< TreeSolver, ReactionNetwork,         // Gillespie for electrolyte
        ReactionNetworkParameters, TrajectoriesSql>
    dispatcher

    LGMC(int argc, char ** argv);
    ~LGMC();
    void print_usage();

};

}

#endif

//#include "mpi.h"
#include "LGMC/LGMC.h"
#include <iostream>

using namespace LGMC_NS;

/* ----------------------------------------------------------------------
   main program to drive RNMCp
------------------------------------------------------------------------- */

int main(int argc, char **argv) {

    LGMC *battery = new LGMC(argc, argv);
    delete battery;
}



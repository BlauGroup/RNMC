/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#ifndef RNMC_SIMULATION_H
#define RNMC_SIMULATION_H

#include <functional>
#include <csignal>
#include <set>
#include <atomic>
#include <unistd.h>
#include <string>
#include <cstring>

#include "../GMC/tree_solver.h"

// In the GNUC Library, sig_atomic_t is a typedef for int,
// which is atomic on all systems that are supported by the
// GNUC Library
volatile sig_atomic_t do_shutdown = 0;

// std::atomic is safe, as long as it is lock-free
std::atomic<bool> shutdown_requested = false;
static_assert(std::atomic<bool>::is_always_lock_free);
// or, at runtime: assert( shutdown_requested.is_lock_free() );

/* ------------------------------------------------------------------- */

template <typename Solver>
class Simulation
{
public:
    unsigned long int seed;
    double time;
    int step; // number of reactions which have occoured
    unsigned long int history_chunk_size;
    std::function<void(Update)> update_function;

    Simulation(unsigned long int seed,
               int history_chunk_size,
               int step,
               double time) : seed(seed),
                              time(time),
                              step(step),
                              history_chunk_size(history_chunk_size) {};

    void execute_steps(int step_cutoff);
    void execute_time(double time_cutoff);
    virtual bool execute_step();
    void write_error_message(std::string s);
};

#include "simulation.cpp"

#endif

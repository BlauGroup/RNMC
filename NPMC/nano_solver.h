/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://lzichi.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#ifndef RNMC_NANO_SOLVER_H
#define RNMC_NANO_SOLVER_H

#include <vector>
#include <optional>
#include <cmath>
#include <map>
#include <csignal>
#include <iostream>

#include "../core/sampler.h"
#include "../core/RNMC_types.h"
#include "NPMC_types.h"

struct NanoUpdate
{
    unsigned long int index;
    NanoReaction reaction; 
};

class NanoSolver
{
private:
    Sampler sampler;
    std::vector<double> cumulative_propensities;
    int number_of_active_indices;
    double propensity_sum;
    unsigned long int last_non_zero_event;

public:
    std::vector<NanoReaction> current_reactions;
    NanoSolver(unsigned long int seed, std::vector<NanoReaction> &current_reactions);
    NanoSolver(unsigned long int seed, std::vector<NanoReaction> &&current_reactions);
    void update();
    void update(NanoUpdate update);
    void update(std::vector<NanoUpdate> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
    NanoSolver() : sampler(Sampler(0)){}; 
};

#endif
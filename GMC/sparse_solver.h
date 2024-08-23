/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#ifndef RNMC_SPARSE_SOLVER_H
#define RNMC_SPARSE_SOLVER_H

#include <vector>
#include <optional>
#include <cmath>
#include <map>

#include "../core/sampler.h"
#include "../core/RNMC_types.h"

class SparseSolver
{
private:
    Sampler sampler;
    std::map<unsigned long int, double> propensities;
    double propensity_sum;

public:
    SparseSolver(unsigned long int seed, std::vector<double> &initial_propensities);
    SparseSolver() : sampler(Sampler(0)){};
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};

#endif
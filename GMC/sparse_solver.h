#ifndef RNMC_SPARSE_SOLVER_H
#define RNMC_SPARSE_SOLVER_H

#include "../core/sampler.h"
#include "../core/RNMC_types.h"

#include <vector>
#include <optional>
#include <cmath>
#include <map>

class SparseSolver {
private:
    Sampler sampler;
    std::map<unsigned long int, double> propensities;
    double propensity_sum;

public:
    SparseSolver(unsigned long int seed, std::vector<double> &initial_propensities);
    SparseSolver(): sampler(Sampler(0)) {}; // defualt constructor
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};

#endif
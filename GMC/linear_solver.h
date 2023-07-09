/* ----------------------------------------------------------------------

    the solver is the algorithmic backbone of a monte carlo simulation
    it decides what will occour next. For now, we have the linear
    solver and a tree solver ported from spparks:
    https://spparks.sandia.gov/

------------------------------------------------------------------------- */

#ifndef RNMC_LINERAR_SOLVER_H
#define RNMC_LINERAR_SOLVER_H

#include "../core/sampler.h"
#include "../core/RNMC_types.h"

#include <vector>
#include <optional>
#include <cmath>
#include <map>

class LinearSolver {
private:
    Sampler sampler;
    std::vector<double> propensities;
    int number_of_active_indices;
    double propensity_sum;
    unsigned long int last_non_zero_event;

public:
    // for linear solver we can moves initial_propensities vector into the object
    // and use it as the propensity buffer. For compatibility with other solvers,
    // we also implement initialization by copying from a reference
    LinearSolver(unsigned long int seed, std::vector<double> &&initial_propensities);
    LinearSolver(unsigned long int seed, std::vector<double> &initial_propensities);
    LinearSolver(): sampler(Sampler(0)) {}; // defualt constructor
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};

#endif
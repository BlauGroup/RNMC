#ifndef RNMC_TREE_SOLVER_H
#define RNMC_TREE_SOLVER_H

#include <vector>
#include <optional>
#include <cmath>
#include <map>

#include "../core/sampler.h"
#include "../core/RNMC_types.h"


// the solver is the algorithmic backbone of a monte carlo simulation
// it decides what will occour next.  for now, we have the linear
// solver and a tree solver ported from spparks:
// https://spparks.sandia.gov/

class TreeSolver {
private:
    Sampler sampler;
    std::vector<double> tree; // we store the propensities in a binary heap
    int number_of_indices; // for this solver, different to length of tree
    int number_of_active_indices; // an index is active if its propensity is non zero
    int propensity_offset; // index where propensities start as leaves of tree

    // walk tree from root to appropriate leaf
    // value is modified when right branch of tree is traversed
    int find_solve_tree(double value);

public:
    // tree solver is constructed using a reference because it ends up
    // forming the tail end of a larger vector
    TreeSolver(): sampler(Sampler(0)) {}; // defualt constructor
    TreeSolver(unsigned long int seed, std::vector<double> &initial_propensities);
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};

#endif
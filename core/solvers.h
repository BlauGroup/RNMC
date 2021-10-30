#pragma once
#include "sampler.h"
#include <vector>
#include <optional>
#include <cmath>

// the solver is the algorithmic backbone of a monte carlo simulation
// it decides what will occour next.  for now, we have the linear
// solver and a tree solver ported from spparks:
// https://spparks.sandia.gov/



struct Update {
    int index;
    double propensity;
};

struct Event {
    int index;
    double dt;
};


class Solver {
public:
    virtual void update(Update update) = 0;
    virtual void update(std::vector<Update> updates) = 0;
    virtual std::optional<Event> event() = 0;
    virtual double get_propensity(int index) = 0;
    virtual double get_propensity_sum() = 0;
};


class LinearSolver : Solver {
private:
    Sampler sampler;
    std::vector<double> propensities;
    int number_of_active_indices;
    double propensity_sum;

public:
    // linear solver moves initial_propensities vector into the object
    // since it can use it unmodified as its propensity buffer
    LinearSolver(unsigned long int seed, std::vector<double> &&initial_propensities);
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};


class TreeSolver : Solver {
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
    TreeSolver(unsigned long int seed, std::vector<double> &initial_propensities);
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};



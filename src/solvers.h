#pragma once
#include "sampler.h"
#include <vector>

struct Update {
    int index;
    double propensity;
};

struct Event {
    int index;
    double dt;
};

class Solver {
    virtual void update(Update update);
    virtual void update(std::vector<Update> updates);
    virtual Event event();
    virtual double get_propensity(int index);
    virtual double get_propensity_sum();
};

class LinearSolver : Solver {
private:
    Sampler sampler;
    std::vector<double> propensities;
    int number_of_active_indices;
    double propensity_sum;

public:
    LinearSolver(unsigned long int seed, std::vector<double> initial_propensities);
    void update(Update update);
    void update(std::vector<Update> updates);
    Event event();
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

public:
    TreeSolver(unsigned long int seed, std::vector<double> initial_propensities);
    void update(Update update);
    void update(std::vector<Update> updates);
    Event event();
    double get_propensity(int index);
    double get_propensity_sum();

};

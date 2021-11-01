#pragma once

#include "reaction_network.h"
#include "../core/solvers.h"

struct HistoryElement {
    int reaction;
    double time;
};

template <typename Solver>
struct Simulation {
    ReactionNetwork &reaction_network;
    unsigned long int seed;
    std::vector<int> state;
    double time;
    int step; // number of reactions which have occoured
    Solver solver;
    std::vector<HistoryElement> history;

    Simulation(ReactionNetwork &reaction_network,
               unsigned long int seed,
               int step_cutoff) :
        reaction_network (reaction_network),
        seed (seed),
        state (reaction_network.initial_state),
        time (0.0),
        step (0),
        solver (seed, reaction_network.initial_propensities),
        history (step_cutoff)
        {};

    bool execute_step();
    void run_for(int step_cutoff);
    bool check_state_positivity();

};

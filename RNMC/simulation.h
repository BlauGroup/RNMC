#pragma once

#include "reaction_network.h"
#include "../core/solvers.h"

struct HistoryElement {
    int reaction;
    double time;
};


struct Simulation {
    ReactionNetwork &reaction_network;
    unsigned long int seed;
    std::vector<int> state;
    double time;
    int step; // number of reactions which have occoured
};

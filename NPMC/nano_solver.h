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

struct NanoParticleParameters {
};

struct Interaction {
    // either 1 site or two site interaction
    int interaction_id;
    int number_of_sites;
    int species_id[2];
    int left_state[2];
    int right_state[2];

    // the units of rate depend on number_of_sites. If number_of_sites = 1, then
    // rate has units 1 / s. If number of sites = 2, then rate has units 1 / s m^6.
    double rate;
};

struct NanoSite {
    double x;
    double y;
    double z;
    int species_id;
};

// For this nano particle simulator, a reaction consists
// of upto two sites and the interaction id.
// if it is an internal interaction, site_id_2 will be -1
// In a reaction, the sites must be within the interaction radius bound.
struct NanoReaction {
    int site_id[2];
    Interaction interaction;

    // rate has units 1 / s
    double rate;
};

struct NanoUpdate {
    unsigned long int index;
    // double propensity;
    NanoReaction reaction;
};

class NanoSolver {
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
    NanoSolver() : sampler(Sampler(0)) {}; // defualt constructor
};



#endif
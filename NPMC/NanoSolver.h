#pragma once
#include "sampler.h"
#include <vector>
#include <optional>
#include <cmath>
#include <map>
#include <csignal>
#include <iostream>

// the solver is the algorithmic backbone of a monte carlo simulation
// it decides what will occour next.  for now, we have the linear
// solver and a tree solver ported from spparks:
// https://spparks.sandia.gov/

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
// For this nano particle simulator, a reaction consists
// of upto two sites and the interaction id.
// if it is an internal interaction, site_id_2 will be -1
// In a reaction, the sites must be within the interaction radius bound.

//Include this here for now, but breakout into own file later
struct Reaction {
    int site_id[2];
    Interaction interaction;

    // rate has units 1 / s
    double rate;
};

struct Update {
    unsigned long int index;
    // double propensity;
    Reaction reaction;
};

struct Event {
    unsigned long int index;
    double dt;
};

class NanoSolver {
private:
    Sampler sampler;
    std::vector<double> cumulative_propensities;
    int number_of_active_indices;
    double propensity_sum;
    unsigned long int last_non_zero_event;

public:
    // for linear solver we can moves initial_propensities vector into the object
    // and use it as the propensity buffer. For compatibility with other solvers,
    // we also implement initialization by copying from a reference
    // LinearSolver(unsigned long int seed, std::vector<double> &&initial_propensities);
    // LinearSolver(unsigned long int seed, std::vector<double> &initial_propensities);
    std::vector<Reaction> current_reactions;
    NanoSolver(unsigned long int seed, std::vector<Reaction> &current_reactions);
    NanoSolver(unsigned long int seed, std::vector<Reaction> &&current_reactions);
    void update();
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};

// LinearSolver implementation
// LinearSolver can operate directly on the passed propensities using a move
NanoSolver::NanoSolver(
    unsigned long int seed,
    std::vector<Reaction> &&current_reactions) :
    sampler (Sampler(seed)),
    cumulative_propensities (current_reactions.size()),
    number_of_active_indices (0),
    propensity_sum (0.0), 
    // if this move isn't here, the semantics is that initial
    // propensities gets moved into a stack variable for the function
    // call and that stack variable is copied into the object.
    current_reactions (std::move(current_reactions)){
        update();
    };
//
NanoSolver::NanoSolver(
    unsigned long int seed,
    std::vector<Reaction> &current_reactions) :
    sampler (Sampler(seed)),
    cumulative_propensities (current_reactions.size()),
    number_of_active_indices (0),
    propensity_sum (0.0),
    current_reactions (current_reactions){
        update();
    };

void NanoSolver::update() {
    cumulative_propensities.resize(current_reactions.size());
    number_of_active_indices = cumulative_propensities.size();
    // std::cerr << "Propensities are: ";
    if (number_of_active_indices > 0) {
        cumulative_propensities[0] = current_reactions[0].rate;
        // std::cerr << current_reactions[0].rate << ", ";
    }

    for (unsigned int i = 1; i < current_reactions.size(); i++) {
      cumulative_propensities[i] = cumulative_propensities[i-1] + current_reactions[i].rate;
    //   std::cerr << current_reactions[i].rate << ", ";
    }
    // std::cerr << "\n";

    propensity_sum = cumulative_propensities[cumulative_propensities.size()-1];
};

// void LinearSolver::update(Update update) {
//
//     if (propensities[update.index] > 0.0) number_of_active_indices--;
//
//     if (update.propensity > 0.0) {
//         number_of_active_indices++;
//         if ( update.index > last_non_zero_event )
//             last_non_zero_event = update.index;
//     }
//
//
//     propensity_sum -= propensities[update.index];
//     propensity_sum += update.propensity;
//     propensities[update.index] = update.propensity;
// };
//
// void LinearSolver::update(std::vector<Update> updates) {
//     for (Update u : updates) {
//         update(u);
//     }
// };

std::optional<Event> NanoSolver::event() {
    if (number_of_active_indices == 0) {
        propensity_sum = 0.0;
        return std::optional<Event>();
    }

    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;

    unsigned long m;

    // This is a linear search. In the future, use a binary search
    for ( m = 0; m < current_reactions.size(); m++){
        if (cumulative_propensities[m] > fraction) break;
    }
    // std::cerr << "Propensity sum: " << propensity_sum << "\n";

    double dt = - std::log(r2) / propensity_sum;
    if (m < current_reactions.size())
        return std::optional<Event> (Event {.index = m, .dt = dt});
    else
        return std::optional<Event> (Event {.index = last_non_zero_event, .dt = dt});
}

double NanoSolver::get_propensity(int index) {
    return current_reactions[index].rate;
}
 
double NanoSolver::get_propensity_sum() {
    return propensity_sum;
}
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

class LinearSolver {
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
    LinearSolver(unsigned long int seed, std::vector<Reaction> &current_reactions);
    LinearSolver(unsigned long int seed, std::vector<Reaction> &&current_reactions);
    void update();
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};


// class TreeSolver {
// private:
//     Sampler sampler;
//     std::vector<double> tree; // we store the propensities in a binary heap
//     int number_of_indices; // for this solver, different to length of tree
//     int number_of_active_indices; // an index is active if its propensity is non zero
//     int propensity_offset; // index where propensities start as leaves of tree
//
//     // walk tree from root to appropriate leaf
//     // value is modified when right branch of tree is traversed
//     int find_solve_tree(double value);
//
// public:
//     // tree solver is constructed using a reference because it ends up
//     // forming the tail end of a larger vector
//     TreeSolver(unsigned long int seed, std::vector<double> &initial_propensities);
//     void update(Update update);
//     void update(std::vector<Update> updates);
//     std::optional<Event> event();
//     double get_propensity(int index);
//     double get_propensity_sum();
// };


// class SparseSolver {
// private:
//     Sampler sampler;
//     std::map<unsigned long int, double> propensities;
//     double propensity_sum;
//
// public:
//     SparseSolver(unsigned long int seed, std::vector<double> &initial_propensities);
//     void update(Update update);
//     void update(std::vector<Update> updates);
//     std::optional<Event> event();
//     double get_propensity(int index);
//     double get_propensity_sum();
// };




// LinearSolver implementation
// LinearSolver can operate directly on the passed propensities using a move
LinearSolver::LinearSolver(
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
LinearSolver::LinearSolver(
    unsigned long int seed,
    std::vector<Reaction> &current_reactions) :
    sampler (Sampler(seed)),
    cumulative_propensities (current_reactions.size()),
    number_of_active_indices (0),
    propensity_sum (0.0),
    current_reactions (current_reactions){
        update();
    };

void LinearSolver::update() {
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

std::optional<Event> LinearSolver::event() {
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

double LinearSolver::get_propensity(int index) {
    return current_reactions[index].rate;
}

double LinearSolver::get_propensity_sum() {
    return propensity_sum;
}


// // TreeSolver implementation
// // TreeSolver always copies the initial propensities into a new array.
// TreeSolver::TreeSolver(
//     unsigned long int seed,
//     std::vector<double> &initial_propensities) :
//     sampler (Sampler(seed)),
//     number_of_active_indices (0) {
//         int m = 0; // tree depth
//         int pow2 = 1; // power of 2 >= numberOfReactions
//         number_of_indices = initial_propensities.size();
//
//         while (pow2 < number_of_indices) {
//             pow2 *= 2;
//             m++;
//         };
//
//         int number_of_tree_nodes = 2 * pow2 - 1;
//         propensity_offset = pow2 - 1;
//         tree.resize(number_of_tree_nodes, 0.0);
//
//         for (int i = propensity_offset;
//              i < propensity_offset + number_of_indices;
//              i++) {
//             tree[i] = initial_propensities[i - propensity_offset];
//             if (tree[i] > 0.0) number_of_active_indices++;
//
//         }
//
//         int child1, child2;
//         for (int parent = propensity_offset - 1; parent >= 0; parent--) {
//             child1 = 2 * parent + 1;
//             child2 = 2 * parent + 2;
//             tree[parent] = tree[child1] + tree[child2];
//         };
//
//
// };
//
// void TreeSolver::update(Update update) {
//     if (tree[propensity_offset + update.index] > 0.0) number_of_active_indices--;
//     if (update.propensity > 0.0) number_of_active_indices++;
//     tree[propensity_offset + update.index] = update.propensity;
//
//     int parent, sibling;
//     int i = propensity_offset + update.index;
//
//     while (i > 0) {
//         if (i % 2) sibling = i + 1;
//         else sibling = i - 1;
//         parent = (i - 1) / 2;
//         tree[parent] = tree[i] + tree[sibling];
//         i = parent;
//     }
// }
//
// void TreeSolver::update(std::vector<Update> updates) {
//     for (Update u : updates)
//         update(u);
// }
//
// int TreeSolver::find_solve_tree(double value) {
//     int i, left_child;
//     i = 0;
//     while (i < propensity_offset) {
//         left_child = 2*i + 1;
//         if (value <= tree[left_child]) i = left_child;
//         else {
//             value -= tree[left_child];
//             i = left_child + 1;
//         }
//     }
//     return i - propensity_offset;
// }
//
// std::optional<Event> TreeSolver::event() {
//     unsigned long int m;
//     double r1,r2, dt;
//
//
//     if (number_of_active_indices == 0) {
//         return std::optional<Event>();
//     }
//
//
//     r1 = sampler.generate();
//     r2 = sampler.generate();
//
//     double value = r1 * tree[0];
//
//     m = find_solve_tree(value);
//     dt = - log(r2) / tree[0];
//
//     return std::optional<Event>(Event {.index = m, .dt = dt});
//
// }
//
// double TreeSolver::get_propensity(int index) {
//     return tree[propensity_offset + index];
// }
//
// double TreeSolver::get_propensity_sum() {
//     return tree[0];
// }


// SparseSolver::SparseSolver(unsigned long int seed, std::vector<double> &initial_propensities) :
//     sampler (Sampler(seed)) {
//
//     for ( int i = 0; i < (int) initial_propensities.size(); i++ ) {
//
//         if ( initial_propensities[i] != 0.0 ) {
//             propensities[i] = initial_propensities[i];
//             propensity_sum += initial_propensities[i];
//         }
//     }
// }
//
// void SparseSolver::update(Update update) {
//     propensity_sum -= propensities[update.index];
//     propensity_sum += update.propensity;
//
//     if ( update.propensity == 0.0 ) {
//         propensities.erase(update.index);
//     } else {
//         propensities[update.index] = update.propensity;
//     }
// }
//
//
// void SparseSolver::update(std::vector<Update> updates) {
//     for ( Update update : updates )
//         this->update(update);
// }
//
//
// std::optional<Event> SparseSolver::event() {
//
//     if (propensities.size() == 0) {
//         return std::optional<Event>();
//     }
//
//     double r1 = sampler.generate();
//     double r2 = sampler.generate();
//     double fraction = propensity_sum * r1;
//     double partial = 0.0;
//
//     std::map<unsigned long int ,double>::iterator it = propensities.begin();
//
//     for ( it = propensities.begin();
//           it != propensities.end();
//           it++
//         ) {
//         partial += std::get<1>(*it);
//         if (partial > fraction) break;
//     }
//
//     double dt = - std::log(r2) / propensity_sum;
//     if ( it != propensities.end() )
//         return std::optional<Event> (Event {.index = std::get<0>(*it), .dt = dt});
//     else
//         return std::optional<Event> (Event {.index = std::get<0>(*propensities.rbegin()), .dt = dt});
// }
//
//
// double SparseSolver::get_propensity(int index) {
//     if ( propensities.find(index) != propensities.end() )
//         return propensities[index];
//     else
//         return 0;
// }
//
//
// double SparseSolver::get_propensity_sum() { return propensity_sum;};

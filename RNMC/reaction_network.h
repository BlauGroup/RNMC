#pragma once
#include <stdint.h>
#include <vector>
#include <optional>
#include <mutex>
#include "../core/sql.h"

struct Reaction {
    // we assume that each reaction has zero, one or two reactants
    uint8_t number_of_reactants;
    uint8_t number_of_products;

    int reactants[2];
    int products[2];

    double rate;
};

struct DependentsNode {

    // reactions which depend on current reaction.
    // dependents will be nothing if it has not been computed.
    std::optional<std::vector<int>> dependents;
    std::mutex mutex;
    int number_of_occurrences; // number of times the reaction has occoured.

    DependentsNode();

};

struct ReactionNetwork {

    std::vector<Reaction> reactions; // list of reactions
    std::vector<int> initial_state; // initial state for all the simulations
    std::vector<double> initial_propensities; // initial propensities for all the reactions

    double factor_zero; // rate modifer for reactions with zero reactants
    double factor_two; // rate modifier for reactions with two reactants
    double factor_duplicate; // rate modifier for reactions of form A + A -> ...

    // number of times a reaction needs to fire before we compute its
    // node in the dependency graph
    int dependency_threshold;

    std::vector<DependentsNode> dependency_graph;

    ReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        int dependency_threshold);

};


#ifndef RNMC_GILLESPIE_REACTION_NETWORK_H
#define RNMC_GILLESPIE_REACTION_NETWORK_H

#include "reaction_network.h"

struct Reaction {
    // we assume that each reaction has zero, one or two reactants
    uint8_t number_of_reactants;
    uint8_t number_of_products;

    int reactants[2];
    int products[2];

    double rate;
    double dG;
};

// parameters passed to the ReactionNetwork constructor
// by the dispatcher which are model specific
struct ReactionNetworkParameters {
};

class GillespieReactionNetwork : public ReactionNetwork {
public:
    uint8_t energy_budget = 0;
    std::vector<Reaction> reactions; // list of reactions

    GillespieReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        ReactionNetworkParameters parameters);

    void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction);
 
};

#endif
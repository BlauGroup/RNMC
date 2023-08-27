#ifndef ENERGY_REACTION_NETWORK_H
#define ENERGY_REACTION_NETWORK_H

#include <stdint.h>
#include <vector>
#include <optional>
#include <mutex>
#include "../core/sql.h"
#include "sql_types.h"
#include "../core/solvers.h"
#include "../core/simulation.h"

struct EnergyReaction {
    // we assume that each reaction has zero, one or two reactants
    uint8_t number_of_reactants;
    uint8_t number_of_products;

    int reactants[2];
    int products[2];

    double rate;
    double dG;
};

// parameters passed to the EnergyReactionNetwork constructor
// by the dispatcher which are model specific
struct EnergyReactionNetworkParameters {
    double energy_budget;
};

class EnergyReactionNetwork : public ReactionNetwork {
public:
    double energy_budget; // The total energy available (for Delta G > 0 reactions)
    std::vector<EnergyReaction> reactions; // list of reactions

    EnergyReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        ReactionNetworkParameters parameters);

    double compute_propensity(
        std::vector<int> &state,
        int reaction,
        double energy_budget);

    void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction);

    void update_energy_budget(
        double &energy_budget,
        int next_reaction);
};



#endif
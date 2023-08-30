#ifndef ENERGY_REACTION_NETWORK_H
#define ENERGY_REACTION_NETWORK_H

#include "reaction_network.h"

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

class EnergyReactionNetwork : public ReactionNetwork<EnergyReaction> {
public:
    double energy_budget; // The total energy available (for Delta G > 0 reactions)

    EnergyReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        EnergyReactionNetworkParameters parameters);

    double compute_energy_propensity(
        std::vector<int> &state,
        int reaction,
        double energy_budget);

    void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction, 
        double energy_budget);

    void update_energy_budget(
        double &energy_budget,
        int next_reaction);
};

/*---------------------------------------------------------------------------*/

EnergyReactionNetwork::EnergyReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     EnergyReactionNetworkParameters parameters)

    {
    
    // call base class constructor 
    ReactionNetwork<EnergyReaction>(reaction_network_database,
                    initial_state_database);

    // Get the energy_budget from parameters
    energy_budget = parameters.energy_budget;

    SqlStatement<EnergyReactionSql> reaction_statement (reaction_network_database);
    SqlReader<EnergyReactionSql> reaction_reader (reaction_statement);

    // reaction_id is lifted so we can do a sanity check after the
    // loop.  Make sure size of reactions vector, last reaction_id and
    // metadata number_of_reactions are all the same

    unsigned long int reaction_id = 0;

    while(std::optional<EnergyReactionSql> maybe_reaction_row = reaction_reader.next()) {

        EnergyReactionSql reaction_row = maybe_reaction_row.value();
        uint8_t number_of_reactants = reaction_row.number_of_reactants;
        uint8_t number_of_products = reaction_row.number_of_products;
        reaction_id = reaction_row.reaction_id;


        EnergyReaction reaction = {
            .number_of_reactants = number_of_reactants,
            .number_of_products = number_of_products,
            .reactants = { reaction_row.reactant_1, reaction_row.reactant_2 },
            .products = { reaction_row.product_1, reaction_row.product_2},
            .rate = reaction_row.rate,
            .dG = reaction_row.dG
        };

        reactions[reaction_id] = reaction;
    }
    
    std::cerr << "energy_budget: " << energy_budget << std::endl;

    std::cerr << time::time_stamp() << "computing dependency graph...\n";
    compute_dependents();
    std::cerr << time::time_stamp() << "finished computing dependency graph\n";

} // EnergyReactionNetwork()

/*---------------------------------------------------------------------------*/

double EnergyReactionNetwork::compute_energy_propensity(
    std::vector<int> &state,
    int reaction_index,
    double energy_budget) {
    // Compute propensities when we are considering dG > 0 reactions

    EnergyReaction &reaction = reactions[reaction_index];

    if (reaction.dG > energy_budget) {
        // When the reaction requires more energy than is available, it cannot happen

        return 0.0;
    } else {
        // If the reaction fits within the energy budget, compute its propensity as usual

        // Note: We allow all dG < 0 reactions to occur as usual.

        return compute_propensity(state,
                                reaction_index);
    }
}

/*---------------------------------------------------------------------------*/

void EnergyReactionNetwork::update_propensities(
    std::function<void(Update update)> update_function,
    std::vector<int> &state,
    int next_reaction, 
    double energy_budget
    ) {

    EnergyReaction &reaction = reactions[next_reaction];

    std::vector<int> species_of_interest;
    species_of_interest.reserve(4);

    for ( int i = 0; i < reaction.number_of_reactants; i++ ) {
        int reactant_id = reaction.reactants[i];
        species_of_interest.push_back(reactant_id);
    }

    for ( int j = 0; j < reaction.number_of_products; j++ ) {
        int product_id = reaction.products[j];
        species_of_interest.push_back(product_id);
    }

    // Update the propensities for reactions corresponding to species which were produced or consumed
    for ( int species_id : species_of_interest ) {
        for ( unsigned int reaction_index : dependents[species_id] ) {
            double new_propensity = compute_energy_propensity(
                state,
                reaction_index,
                energy_budget);

            update_function(Update {
                    .index = reaction_index,
                    .propensity = new_propensity});
        }
    }

    // We'll need to update the propensities for all reactions which have a dG > energy_budget
    for ( unsigned int reaction_index = 0; reaction_index < reactions.size(); reaction_index++){
        double new_propensity = compute_energy_propensity(
            state,
            reaction_index,
            energy_budget);

        update_function(Update {
                .index = reaction_index,
                .propensity = new_propensity});
    }
}

/*---------------------------------------------------------------------------*/

void EnergyReactionNetwork::update_energy_budget(
    double &energy_budget,
    int next_reaction) {

        EnergyReaction &reaction = reactions[next_reaction];

        // We only need to update the energy budget when the reaction triggered is dG > 0.
        // This way we avoid reaction loops
        if ( reaction.dG > 0 ) {
            energy_budget = energy_budget - reaction.dG;
        }
}

#endif
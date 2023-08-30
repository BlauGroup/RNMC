#ifndef RNMC_GILLESPIE_REACTION_NETWORK_H
#define RNMC_GILLESPIE_REACTION_NETWORK_H

#include "reaction_network.h"

struct GillespieReaction {
    // we assume that each reaction has zero, one or two reactants
    uint8_t number_of_reactants;
    uint8_t number_of_products;

    int reactants[2];
    int products[2];

    double rate;
};

class GillespieReactionNetwork : public ReactionNetwork<GillespieReaction> {
public:
    uint8_t energy_budget = 0;

    GillespieReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        ReactionNetworkParameters parameters);

    void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction);
};

#include "gillespie_reaction_network.h"

GillespieReactionNetwork::GillespieReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     ReactionNetworkParameters)
    {

    // call base class constructor 
    ReactionNetwork(reaction_network_database,
    initial_state_database);

    SqlStatement<ReactionSql> reaction_statement (reaction_network_database);
    SqlReader<ReactionSql> reaction_reader (reaction_statement);

    // reaction_id is lifted so we can do a sanity check after the
    // loop.  Make sure size of reactions vector, last reaction_id and
    // metadata number_of_reactions are all the same

    unsigned long int reaction_id = 0;

    while(std::optional<ReactionSql> maybe_reaction_row = reaction_reader.next()) {

        ReactionSql reaction_row = maybe_reaction_row.value();
        uint8_t number_of_reactants = reaction_row.number_of_reactants;
        uint8_t number_of_products = reaction_row.number_of_products;
        reaction_id = reaction_row.reaction_id;

        GillespieReaction reaction = {
            .number_of_reactants = number_of_reactants,
            .number_of_products = number_of_products,
            .reactants = { reaction_row.reactant_1, reaction_row.reactant_2 },
            .products = { reaction_row.product_1, reaction_row.product_2},
            .rate = reaction_row.rate
        };
        reactions[reaction_id] = reaction;
    }
     
    std::cerr << time::time_stamp() << "computing dependency graph...\n";
    compute_dependents();
    std::cerr << time::time_stamp() << "finished computing dependency graph\n";

} // GillespieReactionNetwork()

/*---------------------------------------------------------------------------*/

void GillespieReactionNetwork::update_propensities(
    std::function<void(Update update)> update_function,
    std::vector<int> &state,
    int next_reaction) {

    GillespieReaction &reaction = reactions[next_reaction];

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

    for ( int species_id : species_of_interest ) {
        for ( unsigned int reaction_index : dependents[species_id] ) {

            double new_propensity = compute_propensity(
                state,
                reaction_index);

            update_function(Update {
                    .index = reaction_index,
                    .propensity = new_propensity});
        }
    }
} // update_propensities()


#endif
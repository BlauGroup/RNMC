#pragma once
#include <stdint.h>
#include <vector>
#include <optional>
#include <mutex>
#include "../core/sql.h"
#include "sql_types.h"
#include "../core/solvers.h"
#include "../core/simulation.h"

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
    double energy_budget;
};


struct ReactionNetwork {
    std::vector<Reaction> reactions; // list of reactions
    std::vector<int> initial_state; // initial state for all the simulations
    std::vector<double> initial_propensities; // initial propensities for all the reactions
    double factor_zero; // rate modifer for reactions with zero reactants
    double factor_two; // rate modifier for reactions with two reactants
    double factor_duplicate; // rate modifier for reactions of form A + A -> ...
    double energy_budget; // The total energy available (for Delta G > 0 reactions)

    // maps species to the reactions which involve that species
    std::vector<std::vector<int>> dependents;

    ReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        ReactionNetworkParameters parameters);

    void compute_dependents();

    double compute_propensity(
        std::vector<int> &state,
        int reaction_index
        );

    double compute_propensity(
        std::vector<int> &state,
        int reaction,
        double energy_budget
        );

    void update_state(
        std::vector<int> &state,
        int reaction_index
        );

    std::vector<int> get_species_of_interest(
        Reaction reaction
        );

    void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction
        );
    
    void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction,
        double energy_budget
        );

    void update_energy_budget(
        double &energy_budget,
        int next_reaction
        );
    
    // convert a history element as found a simulation to history
    // to a SQL type.
    TrajectoriesSql history_element_to_sql(
        int seed,
        HistoryElement history_element);

};

ReactionNetwork::ReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     ReactionNetworkParameters parameters)

    {

    // collecting reaction network metadata
    SqlStatement<MetadataSql> metadata_statement (reaction_network_database);
    SqlReader<MetadataSql> metadata_reader (metadata_statement);

    std::optional<MetadataSql> maybe_metadata_row = metadata_reader.next();

    if (! maybe_metadata_row.has_value()) {
        std::cerr << time_stamp()
                  << "no metadata row\n";

        std::abort();
    }

    MetadataSql metadata_row = maybe_metadata_row.value();



    // setting reaction network factors
    SqlStatement<FactorsSql> factors_statement (initial_state_database);
    SqlReader<FactorsSql> factors_reader (factors_statement);

    // TODO: make sure this isn't nothing
    FactorsSql factors_row = factors_reader.next().value();
    factor_zero = factors_row.factor_zero;
    factor_two = factors_row.factor_two;
    factor_duplicate = factors_row.factor_duplicate;

    // Get the energy_budget from parameters
    energy_budget = parameters.energy_budget;

    // loading intial state
    initial_state.resize(metadata_row.number_of_species);

    SqlStatement<InitialStateSql> initial_state_statement (initial_state_database);
    SqlReader<InitialStateSql> initial_state_reader (initial_state_statement);

    int species_id;
    while(std::optional<InitialStateSql> maybe_initial_state_row =
          initial_state_reader.next()) {

        InitialStateSql initial_state_row = maybe_initial_state_row.value();
        species_id = initial_state_row.species_id;
        initial_state[species_id] = initial_state_row.count;
    }

    // loading reactions
    // vectors are default initialized to empty.
    // it is "cleaner" to resize the default vector than to
    // drop it and reinitialize a new vector.
    reactions.resize(metadata_row.number_of_reactions);

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


        Reaction reaction = {
            .number_of_reactants = number_of_reactants,
            .number_of_products = number_of_products,
            .reactants = { reaction_row.reactant_1, reaction_row.reactant_2 },
            .products = { reaction_row.product_1, reaction_row.product_2},
            .rate = reaction_row.rate,
            .dG = reaction_row.dG
        };

        reactions[reaction_id] = reaction;
    
    }

    // sanity check
    if ( metadata_row.number_of_reactions != reaction_id + 1 ||
         metadata_row.number_of_reactions != reactions.size() ) {
        // TODO: improve logging
        std::cerr << time_stamp() <<  "reaction loading failed\n";
        std::abort();
    }

    // computing initial propensities
    initial_propensities.resize(metadata_row.number_of_reactions);
    
    // Trying (and unsuccessfully) using currying to move runtime penalty of energy_budget to compile time
    // // auto compute_propensity_fn = [this](std::vector<int> &state, int reaction_index) -> double { return compute_propensity(state, reaction_index); };
    // if (energy_budget > 0 ) {
    //     auto compute_propensity_fn = [this](std::vector<int> &state, int reaction_index) -> double { return compute_propensity(state, reaction_index, energy_budget); };
    // } else {
    //     auto compute_propensity_fn = [this](std::vector<int> &state, int reaction_index) -> double { return compute_propensity(state, reaction_index); };
    // }
    // for (unsigned long int i = 0; i < initial_propensities.size(); i++) {
    //     initial_propensities[i] = (compute_propensity_fn)(initial_state, i);
    // }

    if ( energy_budget > 0 ) {
        for (unsigned long int i = 0; i < initial_propensities.size(); i++) {
            initial_propensities[i] = compute_propensity(initial_state, i, energy_budget);
        }
    } else {
        for (unsigned long int i = 0; i < initial_propensities.size(); i++) {
            initial_propensities[i] = compute_propensity(initial_state, i);
        }
    }
    

    std::cerr << "energy_budget: " << energy_budget << std::endl;

    std::cerr << time_stamp() << "computing dependency graph...\n";
    compute_dependents();
    std::cerr << time_stamp() << "finished computing dependency graph\n";
};


void ReactionNetwork::compute_dependents() {
    // initializing dependency graph

    dependents.resize(initial_state.size());


    for ( unsigned int reaction_id = 0; reaction_id <  reactions.size(); reaction_id++ ) {
        Reaction &reaction = reactions[reaction_id];

        for ( int i = 0; i < reaction.number_of_reactants; i++ ) {
            int reactant_id = reaction.reactants[i];
            dependents[reactant_id].push_back(reaction_id);
        }


        for ( int j = 0; j < reaction.number_of_products; j++ ) {
            int product_id = reaction.products[j];
            dependents[product_id].push_back(reaction_id);
        }
    }

};


double ReactionNetwork::compute_propensity(
    std::vector<int> &state,
    int reaction_index) {
    // Compute propensities in the standard case

    Reaction &reaction = reactions[reaction_index];

    double p;
    // zero reactants
    if (reaction.number_of_reactants == 0)
        p = factor_zero * reaction.rate;

    // one reactant
    else if (reaction.number_of_reactants == 1)
        p = state[reaction.reactants[0]] * reaction.rate;


    // two reactants
    else {
        if (reaction.reactants[0] == reaction.reactants[1])
            p = factor_duplicate
                * factor_two
                * state[reaction.reactants[0]]
                * (state[reaction.reactants[0]] - 1)
                * reaction.rate;

        else
            p = factor_two
                * state[reaction.reactants[0]]
                * state[reaction.reactants[1]]
                * reaction.rate;
    }

    return p;

};

double ReactionNetwork::compute_propensity(
    std::vector<int> &state,
    int reaction_index,
    double energy_budget) {
    // Compute propensities when we are considering dG > 0 reactions

    Reaction &reaction = reactions[reaction_index];

    if (reaction.dG > energy_budget) {
        // When the reaction requires more energy than is available, it cannot happen

        return 0.0;
    } else {
        // If the reaction fits within the energy budget, compute its propensity as usual

        // Note: We allow all dG < 0 reactions to occur as usual.

        return compute_propensity(state,
                                reaction_index);
    }
};

void ReactionNetwork::update_state(
    std::vector<int> &state,
    int reaction_index) {

    for (int m = 0;
         m < reactions[reaction_index].number_of_reactants;
         m++) {
        state[reactions[reaction_index].reactants[m]]--;
    }

    for (int m = 0;
         m < reactions[reaction_index].number_of_products;
         m++) {
        state[reactions[reaction_index].products[m]]++;
    }

}

std::vector<int> ReactionNetwork::get_species_of_interest(
    Reaction reaction
    ) {
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

    return species_of_interest;
}

void ReactionNetwork::update_propensities(
    std::function<void(Update update)> update_function,
    std::vector<int> &state,
    int next_reaction
    ) {

    Reaction &reaction = reactions[next_reaction];

    std::vector<int> species_of_interest = get_species_of_interest(reaction);

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
}

void ReactionNetwork::update_propensities(
    std::function<void(Update update)> update_function,
    std::vector<int> &state,
    int next_reaction,
    double energy_budget
    ) {

    Reaction &reaction = reactions[next_reaction];

    std::vector<int> species_of_interest = get_species_of_interest(reaction);

    // Update the propensities for reactions corresponding to species which were produced or consumed
    for ( int species_id : species_of_interest ) {
        for ( unsigned int reaction_index : dependents[species_id] ) {
            double new_propensity = compute_propensity(
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
        double new_propensity = compute_propensity(
            state,
            reaction_index,
            energy_budget);

        update_function(Update {
                .index = reaction_index,
                .propensity = new_propensity});
    }
}

void ReactionNetwork::update_energy_budget(
    double &energy_budget,
    int next_reaction
    ) {

        Reaction &reaction = reactions[next_reaction];

        // We only need to update the energy budget when the reaction triggered is dG > 0.
        // This way we avoid reaction loops
        if ( reaction.dG > 0 ) {
            energy_budget = energy_budget - reaction.dG;
        }
        
    }

TrajectoriesSql ReactionNetwork::history_element_to_sql(
    int seed,
    HistoryElement history_element) {
    return TrajectoriesSql {
        .seed = seed,
        .step = history_element.step,
        .reaction_id = history_element.reaction_id,
        .time = history_element.time
    };
}

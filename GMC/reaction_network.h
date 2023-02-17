#pragma once
#include <stdint.h>
#include <vector>
#include <optional>
#include <mutex>
#include "../core/sql.h"
#include "sql_types.h"
#include "../core/solvers.h"
#include "../core/sql_types.h"
#include "../core/simulation.h"

struct Reaction {
    // we assume that each reaction has zero, one or two reactants
    uint8_t number_of_reactants;
    uint8_t number_of_products;

    int reactants[2];
    int products[2];

    double rate;
};

// parameters passed to the ReactionNetwork constructor
// by the dispatcher which are model specific
struct ReactionNetworkParameters {
};


struct ReactionNetwork {
    std::vector<Reaction> reactions; // list of reactions
    std::vector<int> initial_state; // initial state for all the simulations
    std::vector<double> initial_propensities; // initial propensities for all the reactions
    double factor_zero; // rate modifer for reactions with zero reactants
    double factor_two; // rate modifier for reactions with two reactants
    double factor_duplicate; // rate modifier for reactions of form A + A -> ...

    // maps species to the reactions which involve that species
    std::vector<std::vector<int>> dependents;

    ReactionNetwork();

    ReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        ReactionNetworkParameters parameters);

    void compute_dependents();

    double compute_propensity(
        std::vector<int> &state,
        int reaction_index);

    void update_state(
        std::vector<int> &state,
        int reaction_index);

    void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction
        );

    // convert a history element as found a simulation to history
    // to a SQL type.
    ReactionNetworkWriteTrajectoriesSql history_element_to_sql(
        int seed,
        ReactionNetworkTrajectoryHistoryElement history_element);

    ReactionNetworkWriteStateSql state_history_element_to_sql(
        int seed, 
        ReactionNetworkStateHistoryElement history_element);

    WriteCutoffSql cutoff_history_element_to_sql(
        int seed,
        CutoffHistoryElement cutoff_history_element);

    bool read_state(SqlReader<ReactionNetworkReadStateSql> state_reader, 
                    std::map<int, std::vector<int>> &temp_seed_state_map,
                    ReactionNetwork &reaction_network);

    void read_trajectories(SqlReader<ReactionNetworkReadTrajectoriesSql> trajectory_reader, 
                           std::map<int, std::vector<int>> &temp_seed_state_map, 
                           std::map<int, int> &temp_seed_step_map, 
                           std::map<int, double> &temp_seed_time_map,
                           ReactionNetwork &reaction_network);

    void compute_initial_propensities();

    void store_state_history(std::vector<ReactionNetworkStateHistoryElement> &state_packet,
    std::vector<int> &state, ReactionNetwork &reaction_network, unsigned long int &seed);


};

ReactionNetwork::ReactionNetwork() {};

ReactionNetwork::ReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     ReactionNetworkParameters)

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
            .rate = reaction_row.rate
        };

        reactions[reaction_id] = reaction;

        }
     
    initial_propensities.resize(metadata_row.number_of_reactions);

    // sanity check
    if ( metadata_row.number_of_reactions != reaction_id + 1 ||
         metadata_row.number_of_reactions != reactions.size() ) {
        // TODO: improve logging
        std::cerr << time_stamp() <<  "reaction loading failed\n";
        std::abort();
    }

    std::cerr << time_stamp() << "computing dependency graph...\n";
    compute_dependents();
    std::cerr << time_stamp() << "finished computing dependency graph\n";
};

void ReactionNetwork::compute_initial_propensities() {
    // computing initial propensities
    for (unsigned long int i = 0; i < initial_propensities.size(); i++) {
        initial_propensities[i] = compute_propensity(initial_state, i);
    }
}


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


void ReactionNetwork::update_propensities(
    std::function<void(Update update)> update_function,
    std::vector<int> &state,
    int next_reaction
    ) {

    Reaction &reaction = reactions[next_reaction];

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
}


ReactionNetworkWriteTrajectoriesSql ReactionNetwork::history_element_to_sql(
    int seed,
    ReactionNetworkTrajectoryHistoryElement history_element) {
    return ReactionNetworkWriteTrajectoriesSql {
        .seed = seed,
        .step = history_element.step,
        .reaction_id = history_element.reaction_id,
        .time = history_element.time
    };
}

ReactionNetworkWriteStateSql ReactionNetwork::state_history_element_to_sql
    (int seed, ReactionNetworkStateHistoryElement history_element) {
        return ReactionNetworkWriteStateSql {
            .seed = seed,
            .species_id = history_element.species_id,
            .count = history_element.count
        };
}

WriteCutoffSql ReactionNetwork::cutoff_history_element_to_sql(
    int seed,
    CutoffHistoryElement cutoff_history_element) {
        return WriteCutoffSql {
            .seed = seed,
            .step = cutoff_history_element.step,
            .time = cutoff_history_element.time
        };
}

bool ReactionNetwork::read_state(SqlReader<ReactionNetworkReadStateSql> state_reader, 
                std::map<int, std::vector<int>> &temp_seed_state_map,
                ReactionNetwork &reaction_network) {
    
    bool read_interupt_states = false;

    // resize all state vectors 
    for(auto it = temp_seed_state_map.begin(); it != temp_seed_state_map.end(); it++) {
        it->second.resize(reaction_network.initial_state.size());
    }

    while (std::optional<ReactionNetworkReadStateSql> maybe_state_row = state_reader.next()){
        read_interupt_states = true;

        ReactionNetworkReadStateSql state_row = maybe_state_row.value();
        temp_seed_state_map[state_row.seed][state_row.species_id] = state_row.count;
    }

    return read_interupt_states;
}

void ReactionNetwork::read_trajectories(SqlReader<ReactionNetworkReadTrajectoriesSql> trajectory_reader, 
                           std::map<int, std::vector<int>> &temp_seed_state_map, 
                           std::map<int, int> &temp_seed_step_map, 
                           std::map<int, double> &temp_seed_time_map, 
                           ReactionNetwork &reaction_network) {
            
    while (std::optional<ReactionNetworkReadTrajectoriesSql> maybe_trajectory_row = trajectory_reader.next()) {
                        // read_trajectory_states = true;

        ReactionNetworkReadTrajectoriesSql trajectory_row = maybe_trajectory_row.value();
        
        Reaction reaction = reaction_network.reactions[trajectory_row.reaction_id];
        // update reactants
        for (int i = 0; i < reaction.number_of_reactants; i++) {
            temp_seed_state_map[trajectory_row.seed][reaction.reactants[i]] = 
            temp_seed_state_map[trajectory_row.seed][reaction.reactants[i]] - 1;
        }
        // update products
        for (int i = 0; i < reaction.number_of_products; i++) {
            temp_seed_state_map[trajectory_row.seed][reaction.products[i]] = 
            temp_seed_state_map[trajectory_row.seed][reaction.products[i]] + 1;
        }

        if (trajectory_row.step > temp_seed_step_map[trajectory_row.seed]) {
            temp_seed_step_map[trajectory_row.seed] = trajectory_row.step;
            temp_seed_time_map[trajectory_row.seed] = trajectory_row.time;
        }
    }
}

void ReactionNetwork::store_state_history(std::vector<ReactionNetworkStateHistoryElement> &state_packet,
    std::vector<int> &state, ReactionNetwork &reaction_network, unsigned long int &seed) {
    
    for (unsigned int i = 0; i < state.size(); i++) {
                state_packet.push_back(ReactionNetworkStateHistoryElement{
                    .seed = seed,
                    .species_id = static_cast<int>(i),
                    .count = state[i]
                });
    }
}
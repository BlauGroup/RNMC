#ifndef RNMC_REACTION_NETWORK_H
#define RNMC_REACTION_NETWORK_H

#include <stdint.h>
#include <vector>
#include <optional>
#include <mutex>
#include <map>
#include <assert.h>

#include "sql_types.h"
#include "../core/sql.h"
#include "../core/RNMC_types.h"
#include "../core/sql_types.h"
#include "../core/queues.h"

#include <vector>

// parameters passed to the ReactionNetwork constructor
// by the dispatcher which are model specific
struct ReactionNetworkParameters {
};

template <typename Reaction>
class ReactionNetwork {
public:
    std::vector<int> initial_state; // initial state for all the simulations
    std::vector<double> initial_propensities; // initial propensities for all the reactions
    double factor_zero; // rate modifer for reactions with zero reactants
    double factor_two; // rate modifier for reactions with two reactants
    double factor_duplicate; // rate modifier for reactions of form A + A -> ...

    // maps species to the reactions which involve that species
    std::vector<std::vector<int>> dependents;
    std::vector<Reaction> reactions;

    ReactionNetwork();

    ReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database);

    void compute_dependents();

    double compute_propensity(
        std::vector<int> &state,
        int reaction_index);

    void update_state(
        std::vector<int> &state,
        int reaction_index);

    void compute_initial_propensities(std::vector<int> state);

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
    
    void checkpoint(SqlReader<ReactionNetworkReadStateSql> state_reader, 
        SqlReader<ReadCutoffSql> cutoff_reader, 
        SqlReader<ReactionNetworkReadTrajectoriesSql> trajectory_reader, 
        std::map<int, std::vector<int>> &temp_seed_state_map, 
        std::map<int, int> &temp_seed_step_map, 
        SeedQueue &temp_seed_queue, 
        std::map<int, double> &temp_seed_time_map, 
        ReactionNetwork &model);

    void store_checkpoint(std::vector<ReactionNetworkStateHistoryElement> 
        &state_packet, std::vector<int> &state,
        unsigned long int &seed, int step, double time, 
        std::vector<CutoffHistoryElement> &cutoff_packet);
};

template <typename Reaction>
ReactionNetwork<Reaction>::ReactionNetwork() {} // ReactionNetwork()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
ReactionNetwork<Reaction>::ReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database)
    {

    assert(true);

} // ReactionNetwork()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void ReactionNetwork<Reaction>::compute_dependents() {
    // initializing dependency graph
    dependents.resize(initial_state.size());

    for ( unsigned int reaction_id = 0; reaction_id <  reactions.size(); reaction_id++ ) {
        Reaction &reaction = reactions[reaction_id];

        for ( int i = 0; i < reaction.number_of_reactants; i++ ) {
            int reactant_id = reaction.reactants[i];
            if(reaction.number_of_reactants == 1 || (reaction.reactants[0] != reaction.reactants[1])) {
                dependents[reactant_id].push_back(reaction_id);
            }
            else if (i == 0) {
                // if i = 1 then duplicate reactant and don't add dependency twice
                dependents[reactant_id].push_back(reaction_id);
            }
        }
    }
} // compute_dependents()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
double ReactionNetwork<Reaction>::compute_propensity(
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
} //compute_propensity()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void ReactionNetwork<Reaction>::update_state(
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
} // update_state()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void ReactionNetwork<Reaction>::compute_initial_propensities(std::vector<int> state) {
    // computing initial propensities
    for (unsigned long int i = 0; i < initial_propensities.size(); i++) {
        initial_propensities[i] = compute_propensity(state, i);
    }
} // compute_initial_propensities()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
ReactionNetworkWriteTrajectoriesSql ReactionNetwork<Reaction>::history_element_to_sql(
    int seed,
    ReactionNetworkTrajectoryHistoryElement history_element) {
    return ReactionNetworkWriteTrajectoriesSql {
        .seed = seed,
        .step = history_element.step,
        .reaction_id = history_element.reaction_id,
        .time = history_element.time
    };
} //history_element_to_sql()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
ReactionNetworkWriteStateSql ReactionNetwork<Reaction>::state_history_element_to_sql
    (int seed, ReactionNetworkStateHistoryElement history_element) {
        return ReactionNetworkWriteStateSql {
            .seed = seed,
            .species_id = history_element.species_id,
            .count = history_element.count
        };
} // state_history_element_to_sql()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
WriteCutoffSql ReactionNetwork<Reaction>::cutoff_history_element_to_sql(
    int seed, CutoffHistoryElement cutoff_history_element) {
        return WriteCutoffSql {
            .seed = seed,
            .step = cutoff_history_element.step,
            .time = cutoff_history_element.time
        };
} // cutoff_history_element_to_sql()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void ReactionNetwork<Reaction>::checkpoint(SqlReader<ReactionNetworkReadStateSql> state_reader, 
                                    SqlReader<ReadCutoffSql> cutoff_reader, 
                                    SqlReader<ReactionNetworkReadTrajectoriesSql> trajectory_reader, 
                                    std::map<int, std::vector<int>> &temp_seed_state_map, 
                                    std::map<int, int> &temp_seed_step_map, 
                                    SeedQueue &temp_seed_queue, 
                                    std::map<int, double> &temp_seed_time_map, 
                                    ReactionNetwork<Reaction> &model) {
    
    bool read_interrupt_states = false;
    std::vector<int> default_state = model.initial_state;

    while (std::optional<unsigned long int> maybe_seed =
               temp_seed_queue.get_seed()){
                unsigned long int seed = maybe_seed.value();
                temp_seed_state_map.insert(std::make_pair(seed, default_state));
    }

    while (std::optional<ReadCutoffSql> maybe_cutoff_row = cutoff_reader.next()){
            ReadCutoffSql cutoff_row = maybe_cutoff_row.value();

            temp_seed_step_map[cutoff_row.seed] = cutoff_row.step;
            temp_seed_time_map[cutoff_row.seed] = cutoff_row.time;
    }

    while (std::optional<ReactionNetworkReadStateSql> maybe_state_row = state_reader.next()){
        read_interrupt_states = true;

        ReactionNetworkReadStateSql state_row = maybe_state_row.value();
        temp_seed_state_map[state_row.seed][state_row.species_id] = state_row.count;
    }

    if(!read_interrupt_states) {
        while (std::optional<ReactionNetworkReadTrajectoriesSql> maybe_trajectory_row = trajectory_reader.next()) {

            ReactionNetworkReadTrajectoriesSql trajectory_row = maybe_trajectory_row.value();
            
            Reaction reaction = model.reactions[trajectory_row.reaction_id];
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
}

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void ReactionNetwork<Reaction>::store_checkpoint(std::vector<ReactionNetworkStateHistoryElement> 
    &state_packet, std::vector<int> &state,
    unsigned long int &seed, int step, double time, 
    std::vector<CutoffHistoryElement> &cutoff_packet) {
    
    // state information
    for (unsigned int i = 0; i < state.size(); i++) {
        state_packet.push_back(ReactionNetworkStateHistoryElement{
            .seed = seed,
            .species_id = static_cast<int>(i),
            .count = state[i]
        });
    }

    // cutoff information
    cutoff_packet.push_back(CutoffHistoryElement {
        .seed = seed,
        .step = step,
        .time = time
    });
} // store_state_history()

/*---------------------------------------------------------------------------*/

#endif
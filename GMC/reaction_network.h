#ifndef RNMC_REACTION_NETWORK_H
#define RNMC_REACTION_NETWORK_H

#include <stdint.h>
#include <vector>
#include <optional>
#include <mutex>
#include <map>

#include "sql_types.h"
#include "../core/sql.h"
#include "../core/RNMC_types.h"
#include "../core/sql_types.h"
#include "../core/queues.h"

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
        &state_packet, std::vector<int> &state, ReactionNetwork &reaction_network, 
        unsigned long int &seed, int step, double time, 
        std::vector<CutoffHistoryElement> &cutoff_packet);

};

#endif
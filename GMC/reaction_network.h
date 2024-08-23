/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

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

template <typename Reaction>
class ReactionNetwork
{
public:
    double factor_zero;      // rate modifer for reactions with zero reactants
    double factor_two;       // rate modifier for reactions with two reactants
    double factor_duplicate; // rate modifier for reactions of form A + A -> ...

    // maps species to the reactions which involve that species
    std::vector<std::vector<int>> dependents;
    std::vector<Reaction> reactions;

    bool isCheckpoint; // write state, cutoff, trajectories while running or if error

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

    void compute_initial_propensities(std::vector<int> state, std::vector<double> &initial_propensities);

    // convert a history element as found a simulation to history
    // to a SQL type.
    ReactionNetworkWriteTrajectoriesSql history_element_to_sql(
        int seed,
        ReactionNetworkTrajectoryHistoryElement history_element);

    ReactionNetworkWriteStateSql state_history_element_to_sql(
        int seed,
        ReactionNetworkStateHistoryElement history_element);
};

template <typename Reaction>
ReactionNetwork<Reaction>::ReactionNetwork() {} // ReactionNetwork()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void ReactionNetwork<Reaction>::compute_dependents()
{

    for (unsigned int reaction_id = 0; reaction_id < reactions.size(); reaction_id++)
    {
        Reaction &reaction = reactions[reaction_id];

        for (int i = 0; i < reaction.number_of_reactants; i++)
        {
            int reactant_id = reaction.reactants[i];
            if (reaction.number_of_reactants == 1 || (reaction.reactants[0] != reaction.reactants[1]))
            {
                dependents[reactant_id].push_back(reaction_id);
            }
            else if (i == 0)
            {
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
    int reaction_index)
{

    Reaction &reaction = reactions[reaction_index];

    double p;
    // zero reactants
    if (reaction.number_of_reactants == 0)
        p = factor_zero * reaction.rate;

    // one reactant
    else if (reaction.number_of_reactants == 1)
        p = state[reaction.reactants[0]] * reaction.rate;

    // two reactants
    else
    {
        if (reaction.reactants[0] == reaction.reactants[1])
            p = factor_duplicate * factor_two * state[reaction.reactants[0]] * (state[reaction.reactants[0]] - 1) * reaction.rate;

        else
            p = factor_two * state[reaction.reactants[0]] * state[reaction.reactants[1]] * reaction.rate;
    }
    assert(p >= 0);
    return p;
} // compute_propensity()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void ReactionNetwork<Reaction>::update_state(
    std::vector<int> &state,
    int reaction_index)
{

    for (int m = 0;
         m < reactions[reaction_index].number_of_reactants;
         m++)
    {
        state[reactions[reaction_index].reactants[m]]--;
    }

    for (int m = 0;
         m < reactions[reaction_index].number_of_products;
         m++)
    {
        state[reactions[reaction_index].products[m]]++;
    }
} // update_state()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void ReactionNetwork<Reaction>::compute_initial_propensities(std::vector<int> state, 
                                                             std::vector<double> &initial_propensities)
{
    // resize to correct shape
    initial_propensities.resize(reactions.size());

    // computing initial propensities
    for (unsigned long int i = 0; i < initial_propensities.size(); i++)
    {
        initial_propensities[i] = compute_propensity(state, i);
    }

} // compute_initial_propensities()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
ReactionNetworkWriteTrajectoriesSql ReactionNetwork<Reaction>::history_element_to_sql(int seed,
                                                                                      ReactionNetworkTrajectoryHistoryElement history_element)
{
    return ReactionNetworkWriteTrajectoriesSql{
        .seed = seed,
        .step = history_element.step,
        .reaction_id = history_element.reaction_id,
        .time = history_element.time};
} // history_element_to_sql()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
ReactionNetworkWriteStateSql ReactionNetwork<Reaction>::state_history_element_to_sql(int seed, 
                                                                                     ReactionNetworkStateHistoryElement history_element)
{
    return ReactionNetworkWriteStateSql{
        .seed = seed,
        .species_id = history_element.species_id,
        .count = history_element.count};
} // state_history_element_to_sql()

/*---------------------------------------------------------------------------*/

#endif
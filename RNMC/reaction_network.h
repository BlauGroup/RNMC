#pragma once
#include <stdint.h>
#include <vector>
#include <optional>
#include <mutex>
#include "../core/sql.h"
#include "sql_types.h"

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

    DependentsNode() :
    dependents (std::optional<std::vector<int>>()),
    mutex (std::mutex()),
    number_of_occurrences (0) {};
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

    std::optional<std::vector<int>> &get_dependency_node(int reaction_index);
    void compute_dependency_node(int reaction_index);

    double compute_propensity(
        std::vector<int> &state,
        int reaction_index);


};

ReactionNetwork::ReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     int dependency_threshold) :

    dependency_threshold( dependency_threshold ) {

    // collecting reaction network metadata
    SqlStatement<MetadataSql> metadata_statement (reaction_network_database);
    SqlReader<MetadataSql> metadata_reader (metadata_statement);

    // TODO: make sure this isn't nothing
    MetadataSql metadata_row = metadata_reader.next().value();


    // can't resize dependency graph because mutexes are not copyable
    dependency_graph =
        std::vector<DependentsNode> (metadata_row.number_of_reactions);



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

    int reaction_id;

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

    // sanity check
    if ( metadata_row.number_of_reactions != reaction_id + 1 ||
         metadata_row.number_of_reactions != reactions.size() ) {
        // TODO: improve logging
        std::cerr << "reaction loading failed";
        std::abort();
    }

    // computing initial propensities
    initial_propensities.resize(metadata_row.number_of_reactions);
    for (int i = 0; i < initial_propensities.size(); i++) {
        initial_propensities[i] = compute_propensity(initial_state, i);
    }
};


std::optional<std::vector<int>> &ReactionNetwork::get_dependency_node(
    int reaction_index) {

    DependentsNode &node = dependency_graph[reaction_index];

    std::lock_guard(node.mutex);

    if (! node.dependents &&
        node.number_of_occurrences >= dependency_threshold ) {
        compute_dependency_node(reaction_index);
    }

    node.number_of_occurrences++;

    return node.dependents;
};


void ReactionNetwork::compute_dependency_node(int reaction_index) {

    DependentsNode &node = dependency_graph[reaction_index];

    int number_of_dependents_count = 0;
    int j; // reaction index
    int l, m, n; // reactant and product indices

    for (j = 0; j < reactions.size(); j++) {
        bool flag = false;

        for (l = 0; l < reactions[j].number_of_reactants; l++) {
            for (m = 0; m < reactions[reaction_index].number_of_reactants; m++) {
                if (reactions[j].reactants[l] ==
                    reactions[reaction_index].reactants[m])
                    flag = true;
            }

            for (n = 0; n < reactions[reaction_index].number_of_products; n++) {
                if (reactions[j].reactants[l] ==
                    reactions[reaction_index].reactants[n])
                    flag = true;
            }
        }

        if (flag)
            number_of_dependents_count++;
    }

    std::vector<int> dependents (number_of_dependents_count);

    int dependents_counter = 0;
    int current_reaction = 0;
    while (dependents_counter < number_of_dependents_count) {
        bool flag = false;
        for (l = 0;
             l < reactions[current_reaction].number_of_reactants;
             l++) {
            for (m = 0; m < reactions[reaction_index].number_of_reactants; m++) {
                if (reactions[current_reaction].reactants[l] ==
                    reactions[reaction_index].reactants[m])
                    flag = true;
            }

            for (n = 0; n < reactions[reaction_index].number_of_products; n++) {
                if (reactions[current_reaction].reactants[l] ==
                    reactions[reaction_index].products[n])
                    flag = true;
            }
        }

        if (flag) {
            dependents[dependents_counter] = current_reaction;
            dependents_counter++;
        }
        current_reaction++;
    }

    node.dependents = std::optional (std::move(dependents));
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

#include "reaction_network.h"
#include "sql_types.h"

ReactionNetwork::ReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     int dependency_threshold) {

    dependency_threshold = 1;

    // collecting reaction network metadata
    SqlStatement<MetadataSql> metadata_statement (reaction_network_database);
    SqlReader<MetadataSql> metadata_reader (std::ref(metadata_statement));

    // TODO: make sure this isn't nothing
    MetadataSql metadata_row = metadata_reader.next().value();

    // setting reaction network factors
    SqlStatement<FactorsSql> factors_statement (initial_state_database);
    SqlReader<FactorsSql> factors_reader (std::ref(factors_statement));

    // TODO: make sure this isn't nothing
    FactorsSql factors_row = factors_reader.next().value();
    factor_zero = factors_row.factor_zero;
    factor_two = factors_row.factor_two;
    factor_duplicate = factors_row.factor_duplicate;

    // loading intial state
    initial_state.resize(metadata_row.number_of_species);

    SqlStatement<InitialStateSql> initial_state_statement (initial_state_database);
    SqlReader<InitialStateSql> initial_state_reader(
        std::ref(initial_state_statement));

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
    SqlReader<ReactionSql> reaction_reader(std::ref(reaction_statement));


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
};


std::vector<int> &ReactionNetwork::get_dependency_node(
    int reaction_index) {

    DependentsNode &node = std::ref(dependency_graph[reaction_index]);

    std::lock_guard(node.mutex);

    if (! node.dependents &&
        node.number_of_occurrences >= dependency_threshold ) {
        compute_dependency_node(reaction_index);
    }

    node.number_of_occurrences++;

    return std::ref(node.dependents.value());
};


void ReactionNetwork::compute_dependency_node(int reaction_index) {

    DependentsNode &node = std::ref(dependency_graph[reaction_index]);

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

    // TODO: check that this move actually happens as it should
    // the three points in dependents should be zerod out and moved
    // into the dependency graph
    node.dependents = std::optional (std::move(dependents));



};

double ReactionNetwork::compute_propensity(
    std::vector<int> &state,
    int reaction_index) {

    Reaction &reaction = std::ref(reactions[reaction_index]);

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

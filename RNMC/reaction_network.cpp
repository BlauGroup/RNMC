#include "reaction_network.h"
#include "sql_types.h"

ReactionNetwork::ReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     int dependency_threshold) {

    // collecting reaction network metadata
    SqlStatement<MetadataSql> metadata_statement (reaction_network_database);
    SqlReader<MetadataSql> metadata_reader (std::ref(metadata_statement));
    MetadataSql metadata_row = metadata_reader.next().value();

    // vectors are default initialized to empty.
    // it is "cleaner" to resize the default vector than to
    // drop it and reinitialize a new vector.
    reactions.resize(metadata_row.number_of_reactions);

    // setting reaction network factors


    SqlStatement<ReactionSql> reaction_statement (reaction_network_database);
    SqlReader<ReactionSql> reaction_reader(std::ref(reaction_statement));


    // setting reactions attribute
    while(true) {
        // variable is lifted so we can do a sanity check. Make sure
        // size of reactions vector, last reaction_id and metadata
        // number_of_reactions are all the same
        int reaction_id;

        std::optional<ReactionSql> maybe_reaction_row = reaction_reader.next();

        if (maybe_reaction_row) {
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

        else {
            // sanity check
            if ( metadata_row.number_of_reactions != reaction_id + 1 ||
                 metadata_row.number_of_reactions != reactions.size() ) {
                std::cerr << "reaction loading failed";
                std::abort();
            }
            break;
        }
    }

};



DependentsNode &ReactionNetwork::get_dependency_node(
    int reaction_index) {

};


void ReactionNetwork::compute_dependency_node(int reaction_index) {


};

double ReactionNetwork::compute_propensity(
    std::vector<int> &state,
    int reaction_index) {


};

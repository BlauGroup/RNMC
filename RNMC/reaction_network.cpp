#include "reaction_network.h"
#include "sql_types.h"

ReactionNetwork::ReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     int dependency_threshold) {

    SqlReader<MetadataRow> metadata_reader (reaction_network_database);
    MetadataRow metadata_row = metadata_reader.next().value();

    // vectors are default initialized to empty.
    // it is "cleaner" to resize the default vector than to
    // drop it and reinitialize a new vector.
    reactions.resize(metadata_row.number_of_reactions);

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

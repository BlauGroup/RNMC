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

    // SqlConnection sql_connection("./rn.sqlite");
    // SqlReader<ReactionRow> reaction_reader(std::ref(sql_connection));

    // while(true) {
    //     std::optional<ReactionRow> maybe_reaction_row = reaction_reader.next();
    //     if (maybe_reaction_row) {
    //         ReactionRow reaction_row = maybe_reaction_row.value();
    //         std::cout << reaction_row.reaction_id << '\n';
    //     }
    //     else {
    //         break;
    //     }
    // }

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

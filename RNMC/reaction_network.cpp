#include "reaction_network.h"

DependentsNode::DependentsNode() {};

ReactionNetwork::ReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     int dependency_threshold) {


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

#include "reaction_network.h"

int main() {
    SqlConnection reaction_network_database (
        "./test_materials/RNMC/rn.sqlite");

    SqlConnection initial_state_database (
        "./test_materials/RNMC/initial_state.sqlite");

    ReactionNetwork reaction_network (
        std::ref(reaction_network_database),
        std::ref(initial_state_database),
        1);

    std::cout << reaction_network.reactions.size() << '\n';

}

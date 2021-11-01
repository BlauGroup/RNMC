#include "reaction_network.h"
#include "simulation.h"
#include <sqlite3.h>

int main() {
    SqlConnection reaction_network_database (
        "./test_materials/RNMC/rn.sqlite", SQLITE_OPEN_READWRITE);

    SqlConnection initial_state_database (
        "./test_materials/RNMC/initial_state.sqlite", SQLITE_OPEN_READWRITE);

    ReactionNetwork reaction_network (
        std::ref(reaction_network_database),
        std::ref(initial_state_database),
        0);

    std::vector<int> &l = reaction_network.get_dependency_node(0).value();

    for (int i : l) {
        std::cout << i << '\n';
    }

    Simulation<LinearSolver> simulation(std::ref(reaction_network), 42, 200);

    simulation.execute_step();
}

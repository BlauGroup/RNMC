#include "simulation.h"

int main() {
    SqlConnection reaction_network_database (
        "./test_materials/RNMC/rn.sqlite", SQLITE_OPEN_READWRITE);

    SqlConnection initial_state_database (
        "./test_materials/RNMC/initial_state.sqlite", SQLITE_OPEN_READWRITE);

    ReactionNetwork reaction_network (
        reaction_network_database,
        initial_state_database,
        0);

    Simulation<LinearSolver> simulation(reaction_network, 42, 200);

    simulation.execute_steps();
}

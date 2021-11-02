#include "dispatcher.h"

int main() {
    SqlConnection reaction_network_database (
        "./test_materials/RNMC/rn.sqlite", SQLITE_OPEN_READWRITE);

    SqlConnection initial_state_database (
        "./test_materials/RNMC/initial_state.sqlite", SQLITE_OPEN_READWRITE);

    ReactionNetwork reaction_network (
        reaction_network_database,
        initial_state_database,
        0);

    unsigned long int base_seed = 1000;
    unsigned long int number_of_seeds = 1000;
    HistoryQueue history_queue;
    SeedQueue seed_queue (number_of_seeds, base_seed);


    SimulatorPayload<LinearSolver> payload (
        reaction_network,
        history_queue,
        seed_queue,
        200);

    payload.run_simulator();

}

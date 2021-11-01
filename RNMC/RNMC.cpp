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

    Simulation<LinearSolver> simulation (reaction_network, 42, 200);

    simulation.execute_steps(200);

    std::cout << simulation.history.size() << '\n';

    HistoryQueue history_queue;

    history_queue.insert_history(std::move(simulation.history));
    std::vector<HistoryElement> history = history_queue.get_history();
}

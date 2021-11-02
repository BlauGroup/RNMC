#include "dispatcher.h"

int main() {
    std::string reaction_network_database_file = "./test_materials/RNMC/rn.sqlite";
    std::string initial_state_database_file =
        "./test_materials/RNMC/initial_state.sqlite";

    Dispatcher<LinearSolver> dispatcher (
        reaction_network_database_file,
        initial_state_database_file,
        1000,
        1000,
        8,
        200,
        0);

    dispatcher.run_dispatcher();

    std::cout << dispatcher.threads.size();



}

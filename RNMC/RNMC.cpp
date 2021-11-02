#include "dispatcher.h"

int main() {
    std::string reaction_network_database_file = "./test_materials/RNMC/rn.sqlite";
    std::string initial_state_database_file =
        "./test_materials/RNMC/initial_state.sqlite";

    Dispatcher<TreeSolver> dispatcher (
        reaction_network_database_file,
        initial_state_database_file,
        1,
        1000,
        1,
        200,
        1);

    dispatcher.run_dispatcher();




}

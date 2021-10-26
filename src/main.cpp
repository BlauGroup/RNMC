#include "solvers.h"
#include <iostream>

int main() {
    std::vector<double> initial_propensities{0.1, 0.1, 0.2, 0.1, 0.4};
    TreeSolver tree_solver(42, initial_propensities);
    LinearSolver linear_solver(42, std::move(initial_propensities));

    for (int i = 0; i < 10000; i++) {
        Event linear_event = linear_solver.event();
        Event tree_event = tree_solver.event();
        if (linear_event.index != tree_event.index ||
            linear_event.dt != tree_event.dt) {
            std::cout << "non matching event found. index = " << i << '\n';
            std::cout << "linear_event: index = "
                      << linear_event.index
                      << " dt = "
                      << linear_event.dt << '\n';
            std::cout << "tree_event: index = "
                      << tree_event.index
                      << " dt = "
                      << tree_event.dt << '\n';
            break;

        }
    }
}

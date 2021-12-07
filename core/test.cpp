#include "solvers.h"
#include <iostream>


int main() {
    // TODO: make this into a test which passes or fails
    std::vector<double> initial_propensities = {0.1, 0.2, 0.3, 0.1, 0.1};
    TreeSolver tree_solver (42, initial_propensities);
    LinearSolver linear_solver_unsed(42, initial_propensities);
    LinearSolver linear_solver (42, std::move(initial_propensities));

    for (int i = 0; i < 100000; i++) {
        Event linear_event = linear_solver.event().value();
        Event tree_event = tree_solver.event().value();
        if (linear_event.index != tree_event.index) {
            std::cout << "non matching event found." << '\n'
                      << "step = " << i << '\n';

            std::cout << "linear_event: index = "
                      << linear_event.index
                      << ", dt = "
                      << linear_event.dt << '\n';
            std::cout << "tree_event: index = "
                      << tree_event.index
                      << ", dt = "
                      << tree_event.dt << '\n';

            return 1;
        }
    }
    return 0;
}

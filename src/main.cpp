#include "solvers.h"
#include "sql.h"
#include <iostream>


int main() {
    std::vector<double> initial_propensities{0.1, 0.1, 0.2, 0.1, 0.4};
    TreeSolver tree_solver(42, initial_propensities);
    LinearSolver linear_solver(42, std::move(initial_propensities));

    for (int i = 0; i < 100000; i++) {
        Event linear_event = linear_solver.event().value();
        Event tree_event = tree_solver.event().value();
        if (linear_event.index != tree_event.index ||
            linear_event.dt != tree_event.dt) {
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
            break;

        }
    }

    SqlConnection sql_connection("./rn.sqlite");
    SqlReader<ReactionRow> reaction_reader(std::ref(sql_connection));

    while(true) {
        std::optional<ReactionRow> maybe_reaction_row = reaction_reader.next();
        if (maybe_reaction_row) {
            ReactionRow reaction_row = maybe_reaction_row.value();
            std::cout << reaction_row.reaction_id << '\n';
        }
        else {
            break;
        }
    }

}

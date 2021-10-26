#include "solvers.h"
#include <iostream>

int main() {

    std::vector<double> initial_propensities(1000000, 0.1);
    LinearSolver solver(42, std::move(initial_propensities));
    for (double p : initial_propensities) {
        std::cout << p << '\n';
    }

}


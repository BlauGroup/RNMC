#include "solvers.h"
#include <cmath>

LinearSolver::LinearSolver(unsigned long int seed, std::vector<double> initial_propensities) :
    sampler{Sampler(seed)},
    propensities{initial_propensities} {
        int i;
        for (i = 0; i < initial_propensities.size(); i++) {
            propensity_sum += initial_propensities[i];
            if (initial_propensities[i] > 0)
                number_of_active_indices += 1;
        }
};

void LinearSolver::update(Update update) {

    if (propensities[update.index] > 0.0) number_of_active_indices--;
    if (update.propensity > 0.0) number_of_active_indices++;
    propensity_sum -= propensities[update.index];
    propensity_sum += update.propensity;
    propensities[update.index] = update.propensity;
};

void LinearSolver::update(std::vector<Update> updates) {
    for (Update u : updates) {
        update(u);
    }
};

Event LinearSolver::event() {
    if (number_of_active_indices == 0) {
        propensity_sum = 0.0;
        return Event {.index = -1, .dt = 0};
    }

    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;
    double partial = 0.0;

    int m;

    for (m = 0; m < propensities.size(); m++) {
        partial += propensities[m];
        if (partial > fraction) break;
    }

    double dt = - std::log(r2) / propensity_sum;
    if (m < propensities.size())
        return Event {.index = m, .dt = dt};
    else
        return Event {.index = m - 1, .dt = dt};
}

double LinearSolver::get_propensity(int index) {
    return propensities[index];
}

double LinearSolver::get_propensity_sum() {
    return propensity_sum;
}

TreeSolver::TreeSolver(unsigned long int seed, std::vector<double> initial_propensities) :
    sampler{Sampler(seed)} {

        int m = 0; // tree depth
        int pow2 = 1; // power of 2 >= numberOfReactions
        number_of_indices = initial_propensities.size();

        while (pow2 < number_of_indices) {
            pow2 *= 2;
            m++;
        };

        int number_of_tree_nodes = 2 * pow2 - 1;
        propensity_offset = pow2 - 1;
        tree = std::vector<double> (number_of_tree_nodes, 0.0);

        for (int i = propensity_offset; i < propensity_offset + number_of_indices; i++) {
            tree[i] = initial_propensities[i - propensity_offset];
            if (tree[i] > 0.0) number_of_active_indices++;
        }

        propensity_sum = 0.0;

        int child1, child2;
        for (int parent = propensity_offset - 1; parent >= 0; parent--) {
            child1 = 2 * parent + 1;
            child2 = 2 * parent + 2;
            tree[parent] = tree[child1] + tree[child2];
        };

        propensity_sum = tree[0];

};


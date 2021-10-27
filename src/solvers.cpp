#include "solvers.h"
#include <cmath>


// LinearSolver implementation
// LinearSolver can opperate directly on the passed propensities using a move
LinearSolver::LinearSolver(
    unsigned long int seed,
    std::vector<double> initial_propensities) :
    sampler{Sampler(seed)},
    propensities{initial_propensities},
    number_of_active_indices{0},
    propensity_sum{0.0} {
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
        return Event {.index = -1, .dt = 0, .event_occoured = false};
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
        return Event {.index = m, .dt = dt, .event_occoured = true};
    else
        return Event {.index = m - 1, .dt = dt, .event_occoured = true};
}

double LinearSolver::get_propensity(int index) {
    return propensities[index];
}

double LinearSolver::get_propensity_sum() {
    return propensity_sum;
}


// TreeSolver implementation
// TreeSolver always copies the initial propensities into a new array.
TreeSolver::TreeSolver(
    unsigned long int seed,
    std::vector<double> initial_propensities) :
    sampler{Sampler(seed)},
    number_of_active_indices{0} {
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

        for (int i = propensity_offset;
             i < propensity_offset + number_of_indices;
             i++) {
            tree[i] = initial_propensities[i - propensity_offset];
            if (tree[i] > 0.0) number_of_active_indices++;

        }

        int child1, child2;
        for (int parent = propensity_offset - 1; parent >= 0; parent--) {
            child1 = 2 * parent + 1;
            child2 = 2 * parent + 2;
            tree[parent] = tree[child1] + tree[child2];
        };


};

void TreeSolver::update(Update update) {
    if (tree[propensity_offset + update.index] > 0.0) number_of_active_indices--;
    if (update.propensity > 0.0) number_of_active_indices++;
    tree[propensity_offset + update.index] = update.propensity;

    int parent, sibling;
    int i = propensity_offset + update.index;

    while (i > 0) {
        if (i % 2) sibling = i + 1;
        else sibling = i - 1;
        parent = (i - 1) / 2;
        tree[parent] = tree[i] + tree[sibling];
        i = parent;
    }
}

void TreeSolver::update(std::vector<Update> updates) {
    for (Update u : updates)
        update(u);
}

int TreeSolver::find_solve_tree(double value) {
    int i, left_child;
    i = 0;
    while (i < propensity_offset) {
        left_child = 2*i + 1;
        if (value <= tree[left_child]) i = left_child;
        else {
            value -= tree[left_child];
            i = left_child + 1;
        }
    }
    return i - propensity_offset;
}

Event TreeSolver::event() {
    int m;
    double r1,r2, dt;


    if (number_of_active_indices == 0) {
        return Event {.index = -1, .dt = 0.0, .event_occoured = false};
    }


    r1 = sampler.generate();
    r2 = sampler.generate();

    double value = r1 * tree[0];

    m = find_solve_tree(value);
    dt = - log(r2) / tree[0];

    return Event {.index = m, .dt = dt, .event_occoured = true};

}

double TreeSolver::get_propensity(int index) {
    return tree[propensity_offset + index];
}

double TreeSolver::get_propensity_sum() {
    return tree[0];
}

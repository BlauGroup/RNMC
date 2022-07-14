#pragma once
#include "sampler.h"
#include <vector>
#include <optional>
#include <cmath>
#include <map>

// the solver is the algorithmic backbone of a monte carlo simulation
// it decides what will occour next.  for now, we have the linear
// solver and a tree solver ported from spparks:
// https://spparks.sandia.gov/

struct Update {
    unsigned long int index;
    double propensity;
};

struct Event {
    unsigned long int index;
    double dt;
};

class LinearSolver {
private:
    Sampler sampler;
    std::vector<double> propensities;
    int number_of_active_indices;
    double propensity_sum;
    unsigned long int last_non_zero_event;

public:
    // for linear solver we can moves initial_propensities vector into the object
    // and use it as the propensity buffer. For compatibility with other solvers,
    // we also implement initialization by copying from a reference
    template < typename Model>
    LinearSolver(unsigned long int seed, Model &&model);

    template < typename Model>
    LinearSolver(unsigned long int seed, Model &model);
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};


class TreeSolver {
private:
    Sampler sampler;
    std::vector<double> tree; // we store the propensities in a binary heap
    int number_of_indices; // for this solver, different to length of tree
    int number_of_active_indices; // an index is active if its propensity is non zero
    int propensity_offset; // index where propensities start as leaves of tree

    // walk tree from root to appropriate leaf
    // value is modified when right branch of tree is traversed
    int find_solve_tree(double value);

public:
    // tree solver is constructed using a reference because it ends up
    // forming the tail end of a larger vector
    template < typename Model>
    TreeSolver(unsigned long int seed, Model &model);
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};


class SparseSolver {
private:
    Sampler sampler;
    std::map<unsigned long int, double> propensities;
    double propensity_sum;

public:
    template < typename Model>
    SparseSolver(unsigned long int seed, Model &model);
    void update(Update update);
    void update(std::vector<Update> updates);
    std::optional<Event> event();
    double get_propensity(int index);
    double get_propensity_sum();
};

// LinearSolver implementation
// LinearSolver can opperate directly on the passed propensities using a move
template < typename Model>
LinearSolver::LinearSolver(
    unsigned long int seed,
    Model &&model) :
    sampler (Sampler(seed)),
    // if this move isn't here, the semantics is that initial
    // propensities gets moved into a stack variable for the function
    // call and that stack variable is copied into the object.
    propensities (std::move(model.initial_propensities)),
    number_of_active_indices (0),
    propensity_sum (0.0) {
        for (unsigned long i = 0; i < propensities.size(); i++) {
            propensity_sum += propensities[i];
            if (propensities[i] > 0) {
                number_of_active_indices += 1;
                last_non_zero_event = i;
            }

        }
    };

template < typename Model>
LinearSolver::LinearSolver(
    unsigned long int seed,
    Model &model) :
    sampler (Sampler(seed)),
    propensities (model.initial_propensities),
    number_of_active_indices (0),
    propensity_sum (0.0) {

        for (unsigned long i = 0; i < propensities.size(); i++) {
            propensity_sum += propensities[i];
            if (propensities[i] > 0) {
                number_of_active_indices += 1;
                last_non_zero_event = i;
            }
        }
    };


void LinearSolver::update(Update update) {

    if (propensities[update.index] > 0.0) number_of_active_indices--;

    if (update.propensity > 0.0) {
        number_of_active_indices++;
        if ( update.index > last_non_zero_event )
            last_non_zero_event = update.index;
    }


    propensity_sum -= propensities[update.index];
    propensity_sum += update.propensity;
    propensities[update.index] = update.propensity;
};

void LinearSolver::update(std::vector<Update> updates) {
    for (Update u : updates) {
        update(u);
    }
};

std::optional<Event> LinearSolver::event() {
    if (number_of_active_indices == 0) {
        propensity_sum = 0.0;
        return std::optional<Event>();
    }

    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;
    double partial = 0.0;

    unsigned long m;

    for (m = 0; m < propensities.size(); m++) {
        partial += propensities[m];
        if (partial > fraction) break;
    }

    double dt = - std::log(r2) / propensity_sum;
    if (m < propensities.size())
        return std::optional<Event> (Event {.index = m, .dt = dt});
    else
        return std::optional<Event> (Event {.index = last_non_zero_event, .dt = dt});
}

double LinearSolver::get_propensity(int index) {
    return propensities[index];
}

double LinearSolver::get_propensity_sum() {
    return propensity_sum;
}


// TreeSolver implementation
// TreeSolver always copies the initial propensities into a new array.
template <typename Model>
TreeSolver::TreeSolver(
    unsigned long int seed,
    Model &model) :
    sampler (Sampler(seed)),
    number_of_active_indices (0) {
        int m = 0; // tree depth
        int pow2 = 1; // power of 2 >= numberOfReactions
        number_of_indices = model.initial_propensities.size();

        while (pow2 < number_of_indices) {
            pow2 *= 2;
            m++;
        };

        int number_of_tree_nodes = 2 * pow2 - 1;
        propensity_offset = pow2 - 1;
        tree.resize(number_of_tree_nodes, 0.0);

        for (int i = propensity_offset;
             i < propensity_offset + number_of_indices;
             i++) {
            tree[i] = model.initial_propensities[i - propensity_offset];
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

std::optional<Event> TreeSolver::event() {
    unsigned long int m;
    double r1,r2, dt;


    if (number_of_active_indices == 0) {
        return std::optional<Event>();
    }


    r1 = sampler.generate();
    r2 = sampler.generate();

    double value = r1 * tree[0];

    m = find_solve_tree(value);
    dt = - log(r2) / tree[0];

    return std::optional<Event>(Event {.index = m, .dt = dt});

}

double TreeSolver::get_propensity(int index) {
    return tree[propensity_offset + index];
}

double TreeSolver::get_propensity_sum() {
    return tree[0];
}

template < typename Model>
SparseSolver::SparseSolver(unsigned long int seed, Model &model) :
    sampler (Sampler(seed)) {

    for ( int i = 0; i < (int) model.initial_propensities.size(); i++ ) {

        if ( model.initial_propensities[i] != 0.0 ) {
            propensities[i] = model.initial_propensities[i];
            propensity_sum += model.initial_propensities[i];
        }
    }
}

void SparseSolver::update(Update update) {
    propensity_sum -= propensities[update.index];
    propensity_sum += update.propensity;

    if ( update.propensity == 0.0 ) {
        propensities.erase(update.index);
    } else {
        propensities[update.index] = update.propensity;
    }
}


void SparseSolver::update(std::vector<Update> updates) {
    for ( Update update : updates )
        this->update(update);
}


std::optional<Event> SparseSolver::event() {

    if (propensities.size() == 0) {
        return std::optional<Event>();
    }

    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;
    double partial = 0.0;

    std::map<unsigned long int ,double>::iterator it = propensities.begin();

    for ( it = propensities.begin();
          it != propensities.end();
          it++
        ) {
        partial += std::get<1>(*it);
        if (partial > fraction) break;
    }

    double dt = - std::log(r2) / propensity_sum;
    if ( it != propensities.end() )
        return std::optional<Event> (Event {.index = std::get<0>(*it), .dt = dt});
    else
        return std::optional<Event> (Event {.index = std::get<0>(*propensities.rbegin()), .dt = dt});
}


double SparseSolver::get_propensity(int index) {
    if ( propensities.find(index) != propensities.end() )
        return propensities[index];
    else
        return 0;
}


double SparseSolver::get_propensity_sum() { return propensity_sum;};




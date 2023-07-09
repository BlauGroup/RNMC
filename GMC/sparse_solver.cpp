#include "sparse_solver.h"

SparseSolver::SparseSolver(unsigned long int seed, 
                            std::vector<double> &initial_propensities) :
    sampler (Sampler(seed)) {

    for ( int i = 0; i < (int) initial_propensities.size(); i++ ) {

        if ( initial_propensities[i] != 0.0 ) {
            propensities[i] = initial_propensities[i];
            propensity_sum += initial_propensities[i];
        }
    }
}

/*---------------------------------------------------------------------------*/

void SparseSolver::update(Update update) {
    propensity_sum -= propensities[update.index];
    propensity_sum += update.propensity;

    if ( update.propensity == 0.0 ) {
        propensities.erase(update.index);
    } else {
        propensities[update.index] = update.propensity;
    }
}

/*---------------------------------------------------------------------------*/

void SparseSolver::update(std::vector<Update> updates) {
    for ( Update update : updates )
        this->update(update);
}

/*---------------------------------------------------------------------------*/

std::optional<Event> SparseSolver::event() {

    if (propensities.size() == 0) {
        return std::optional<Event>();
    }

    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;
    double partial = 0.0;

    std::map<unsigned long int, double>::iterator it = propensities.begin();

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

/*---------------------------------------------------------------------------*/

double SparseSolver::get_propensity(int index) {
    if ( propensities.find(index) != propensities.end() )
        return propensities[index];
    else
        return 0;
}

/*---------------------------------------------------------------------------*/

double SparseSolver::get_propensity_sum() { return propensity_sum;}

/*---------------------------------------------------------------------------*/
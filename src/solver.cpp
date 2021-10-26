#include "solvers.h"
#include <cmath>

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

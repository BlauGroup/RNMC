#ifndef LAT_SOLVER_H
#define LAT_SOLVER_H

#include "../core/sampler.h"
#include <vector>
#include <unordered_map>
#include "../core/solvers.h"
#include <string>
#include <numeric>

class LatSolver {
    private:

        int number_of_active_indices;               // end simulation of no sites with non zero propensity

        double propensity_sum;              

        unsigned long int last_non_zero_event;   

        Sampler sampler; `                          // random number generator

        Model *model;                               // pointer to model 

    public:

        template < typename Model>

        LatSolver(unsigned long int seed, Model &model);

        void update(Update update);

        void update(std::vector<Update> updates);

        std::optional<Event> event();

};



/* ---------------------------------------------------------------------- */

template < typename Model>
LatSolver::LatSolver(unsigned long int seed, Model &model) {
    
    sampler = Sampler(seed);
    number_of_active_indices = 0;
    model = model;

    // lattice propensities initially all zero
    propensity_sum = std::accumulate(model.react_net->propensities.begin(), 
                     model.react_net->propensities.end(), 0);
}


void LatSolver::update(Update update) {
    
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

void LatSolver::update(std::vector<Update> updates) {
    for (Update u : updates) {
        update(u);
    }
};

std::optional<Event> LatSolver::event() {
    if (number_of_active_indices == 0) {
        propensity_sum = 0.0;
        return std::optional<Event>();
    }

    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;
    double partial = 0.0;

    unsigned long m;
    bool isFound = false;
    bool isLattice = false;
    int reaction_id = 0;
    int site_one = std::optional<int>();
    int site_two = 0;
    std::string hash;

    // model handles updating lattice propensities
    propensity_sum -= model.prop_sum;

    // start with Gillespie propensities
    for (m = 0; m < model->react_net->propensities.size(); m++) {
        partial += model->react_net->propensities[m];
        if (partial > fraction) {
            isFound = true;
            reaction_id = m;
            break;
        }
    }

    // go through lattice propensities if not found 
    if(!isFound) {
        for(auto it = model.props.begin(); it != model.props.end(); it++) {
            partial += it->second->first;

            if(partial > fraction) {
                isFound = true;
                isLattice = true;
                hash = it->first;
                reaction_id = it->second->second;

                std::size_t pos = hash.find(".");
                site_one = std::optional<int>(int(hash.substr(0, pos)));
                site_two = std::optional<int>(int(hash.substr(pos)));
                break;
            }

        }
    }
    
    if(isFound ) {
        model.update(site_one, site_two, reaction_id);
    }
    double dt = - std::log(r2) / propensity_sum;

    // model handles updating lattice propensities
    propensity_sum += model.prop_sum;

    if (isFound)
        return std::optional<Event> (Event {.index = m, .dt = dt});
    else
        return std::optional<Event> (Event {.index = last_non_zero_event, .dt = dt});
}




#endif
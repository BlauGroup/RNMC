#ifndef LAT_SOLVER_H
#define LAT_SOLVER_H

#include "../core/sampler.h"
#include <vector>
#include <unordered_map>
#include "../core/solvers.h"
#include <string>
#include <numeric>
#include "LGMC.h"

using namespace LGMC_NS;


struct LatticeEvent {
    int site_one;
    int site_two;
    unsigned long int index;
    double dt;
};

struct LatticeUpdate {
    int site_one;
    int site_two;
    double propensity;
};

class LatSolver {
    public:
        LatSolver(unsigned long int seed, std::vector<double> &&initial_propensities);
        LatSolver(unsigned long int seed, std::vector<double> &initial_propensities);

        void update(Update update);

        void update(LatticeUpdate lattice_update);

        void update(std::vector<Update> updates);

        update(LatticeUpdate lattice_update);

        int sum_row(std::string hash);

        std::pair<std::optional<Event>, std::optional<LatticeEvent>> event();

        void init_lattice_props();

        std::string make_string(int site_one, int site_two);

    private:

        int number_of_active_indices;               // end simulation of no sites with non zero propensity

        double propensity_sum;              

        unsigned long int last_non_zero_event;   

        std::unordered_map<std::string,                     // lattice propensities 
        std::vector< std::pair<double, int> > > props;      // key: (site_one, site_two) value: propensity 
        std::vector<double> propensities;                   // Gillepsie propensities 
};

/* ---------------------------------------------------------------------- */

// LinearSolver implementation
// LinearSolver can opperate directly on the passed propensities using a move
LatSolver::LatSolver(unsigned long int seed,
    std::vector<double> &&initial_propensities) :
    sampler (Sampler(seed)),
    // if this move isn't here, the semantics is that initial
    // propensities gets moved into a stack variable for the function
    // call and that stack variable is copied into the object.
    propensities (std::move(initial_propensities)),
    number_of_active_indices (0),
    propensity_sum (0.0) 
{
    for (unsigned long i = 0; i < propensities.size(); i++) {
            propensity_sum += propensities[i];
            if (propensities[i] > 0) {
                number_of_active_indices += 1;
                last_non_zero_event = i;
            }

        }
};

/* ---------------------------------------------------------------------- */

LatSolver::LatSolver( unsigned long int seed,
    std::vector<double> &initial_propensities) :
    sampler (Sampler(seed)),
    propensities (initial_propensities),
    number_of_active_indices (0),
    propensity_sum (0.0) 
    {
        for (unsigned long i = 0; i < propensities.size(); i++) {
            propensity_sum += propensities[i];
            if (propensities[i] > 0) {
                number_of_active_indices += 1;
                last_non_zero_event = i;
            }

        }
    };

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

void LatSolver::update(LatticeUpdate lattice_update) {

    std::string hash = make_string(lattice_update.site_one, lattice_update.site_two);

    // add propensity 
    props[site_combo].push_back(std::make_pair(lattice_update.propensity, react_id));

    // update running sum 
    propensity_sum += lattice_update.propensity;


    if (update.propensity > 0.0) {
        number_of_active_indices++;
        if ( update.index > last_non_zero_event )
            last_non_zero_event = update.index;
    }

    // hash should only store propensities > 0
    number_of_active_indices++;

};

/* ---------------------------------------------------------------------- */

void LatSolver::update(std::vector<LatticeUpdate> lattice_updates) {
    for (Update u : lattice_updates) {
        update(u);
    }
};

/* ---------------------------------------------------------------------- */

void LatSolver::update(std::vector<Update> updates) {
    for (Update u : updates) {
        update(u);
    }
};

/* ---------------------------------------------------------------------- */

std::pair<std::optional<Event>, std::optional<LatticeEvent>> LatSolver::event() {
    if (number_of_active_indices == 0) {
        propensity_sum = 0.0;
        return std::pair<std::optional<Event>(), std::optional<LatticeEvent>()>;
    }

    double r1 = lgmc->sampler.generate();
    double r2 = lgmc->sampler.generate();
    double fraction = propensity_sum * r1;
    double partial = 0.0;

    unsigned long m;
    bool isFound = false;
    bool isLattice = false;
    int reaction_id = 0;
    std::optional<int> site_one = std::optional<int>();
    std::optional<int> site_two = std::optional<int>();
    std::string hash;

    // model handles updating lattice propensities
    propensity_sum -= lgmc->prop_sum;

    // start with Gillespie propensities
    for (m = 0; m < lgmc->react_net->initial_propensities.size(); m++) {
        partial += lgmc->react_net->initial_propensities[m];
        if (partial > fraction) {
            isFound = true;
            reaction_id = m;
            break;
        }
    }

    // go through lattice propensities if not found 
    if(!isFound) {
        auto it = lgmc->props.begin();
        while(!isFound || it != lgmc->props.end()) {
            for(int i = 0; i < it->second.size(); i++ ) {
                partial += it->second[i].first;

                if(partial > fraction) {
                    isFound = true;
                    isLattice = true;
                    hash = it->first;
                    reaction_id = it->second[i].second;

                    std::size_t pos = hash.find(".");
                    site_one = std::optional<int>(stoi(hash.substr(0, pos)));
                    site_two = std::optional<int>(stoi(hash.substr(pos)));
                    isFound = true;
                    break;
                }
            }
            it++;
        }
    }
    
    double dt = - std::log(r2) / propensity_sum;

    if(isFound && isLattice) {
        return std::pair<std::optional<Event>(), std::optional<LatticeEvent> ( LatticeEvent {.index = reaction_id,
                                                                                .site_one = site_one, .site_two = site_two, 
                                                                                .dt = dt})>;
    }
    else if(isFound && !isLattice) {
        return std::pair<std::optional<Event>(Event {.index = last_non_zero_event, .dt = dt}), std::optional<LatticeEvent> ()>;
    }
    else {
        return std::pair<std::optional<Event> (Event {.index = last_non_zero_event, .dt = dt}), std::optional<LatticeEvent> ()>;
    }
        
}

/* ---------------------------------------------------------------------- */

int LatSolver::sum_row(std::string hash) {
    
    int sum = 0;    
    for(auto it = props[hash].begin(); it != props[hash].end(); it++) {
        sum += it->first;
    }
    return sum;
}

/* ---------------------------------------------------------------------- */

std::string LGMC::make_string(int site_one, int site_two) {

    return (site_one < site_two) ? std::to_string(site_one) + "." + 
    std::to_string(site_two) : std::to_string(site_two) + "." + std::to_string(site_one);

} // make_string

/* ---------------------------------------------------------------------- */

void LatSolver::init_lattice_props() {

} // init_lattice_props()

#endif
#pragma once
#include "../core/sampler.h"
#include <vector>
#include <unordered_map>
#include "../core/solvers.h"
#include <string>
#include <numeric>

struct LatticeUpdate {
    unsigned long int index;
    double propensity;
    int site_one;
    int site_two;
};

struct LatticeEvent {
    std::optional<int> site_one;
    std::optional<int> site_two;
    unsigned long int index;
    double dt;
};

class LatSolver {
    public:
        LatSolver() : sampler(Sampler(0)) {}; // defualt constructor
        LatSolver(unsigned long int seed, std::vector<double> &&initial_propensities);
        LatSolver(unsigned long int seed, std::vector<double> &initial_propensities);

        void update(Update update);
        void update(std::vector<Update> updates);

        void update(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props);
        void update(std::vector<LatticeUpdate> lattice_updates, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props);

        std::optional<LatticeEvent> event_lattice(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props);

        std::string make_string(int site_one, int site_two);

        // for comptability 
        std::optional<Event> event() {return std::optional<Event>();};

        long double propensity_sum;
        int number_of_active_indices;               // end simulation of no sites with non zero propensity   

    private:
        Sampler sampler;
        std::vector<double> propensities;                   // Gillepsie propensities            

};

/* ---------------------------------------------------------------------- */

// LinearSolver implementation
// LinearSolver can opperate directly on the passed propensities using a move
LatSolver::LatSolver(unsigned long int seed,
    std::vector<double> &&initial_propensities) :
    propensity_sum (0.0),
    sampler (Sampler(seed)),
    // if this move isn't here, the semantics is that initial
    // propensities gets moved into a stack variable for the function
    // call and that stack variable is copied into the object.
    propensities (std::move(initial_propensities)),
    number_of_active_indices (0)
{
    for (unsigned long i = 0; i < propensities.size(); i++) {
            propensity_sum += propensities[i];
            if (propensities[i] > 0) {
                number_of_active_indices += 1;
            }

        }
};

/* ---------------------------------------------------------------------- */

LatSolver::LatSolver( unsigned long int seed,
    std::vector<double> &initial_propensities) :
    propensity_sum (0.0),
    sampler (Sampler(seed)),
    propensities (initial_propensities),
    number_of_active_indices (0)
    {
        for (unsigned long i = 0; i < propensities.size(); i++) {
            propensity_sum += propensities[i];
            if (propensities[i] > 0) {
                number_of_active_indices += 1;
            }

        }
    };

/* ---------------------------------------------------------------------- */

void LatSolver::update(Update update) {

    if (propensities[update.index] > 0.0) number_of_active_indices--;

    if (update.propensity > 0.0) {
        number_of_active_indices++;
    }


    propensity_sum -= propensities[update.index];
    propensity_sum += update.propensity;
    propensities[update.index] = update.propensity;

    if(propensity_sum < 0) {
        assert(false);
    }

};

/* ---------------------------------------------------------------------- */

void LatSolver::update(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props) {
    
    propensity_sum += lattice_update.propensity;
    number_of_active_indices++;

    std::string hash = make_string(lattice_update.site_one, lattice_update.site_two);
    props[hash].push_back(std::make_pair(lattice_update.propensity, lattice_update.index));

    double sum = 0;
    for(auto it : props){
        for(int i =0; i < it.second.size(); i++){
            sum += it.second[i].second;
        }
    }
    for(int i = 0; i < propensities.size(); i++){
        sum += propensities[i];
    }

    if(sum != propensity_sum){
        std::cout << "ERROR" << std::endl;
        propensity_sum = sum;
    }

};

/* ---------------------------------------------------------------------- */

void LatSolver::update(std::vector<LatticeUpdate> lattice_updates, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props) {
    for (LatticeUpdate u : lattice_updates) {
        update(u, props);
    }
};

/* ---------------------------------------------------------------------- */

void LatSolver::update(std::vector<Update> updates) {
    for (Update u : updates) {
        update(u);
    }
};

/* ---------------------------------------------------------------------- */

std::optional<LatticeEvent> LatSolver::event_lattice(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props) {
    if (number_of_active_indices == 0) {
        propensity_sum == 0.0;
        return std::optional<LatticeEvent> ();
    }
    
    bool isFound = false;
    unsigned long m;
    unsigned long int reaction_id = 0;
    std::optional<int> site_one = std::optional<int>();
    std::optional<int> site_two = std::optional<int>();
    std::string hash;
    double dt = 0.;


    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;
    long double partial = 0.0;

    // start with Gillespie propensities
    for (m = 0; m < propensities.size(); m++) {
        partial += propensities[m];
        if (partial > fraction) {
            isFound = true;
            reaction_id = m;
            break;
        }
    }

    // go through lattice propensities if not found 
    if(!isFound) {
        auto it = props.begin();
        while(!isFound && it != props.end()) {
            for(int i = 0; i < int(it->second.size()); i++ ) {
                
                partial += it->second[i].first;

                if(partial > fraction) {

                    isFound = true;
                    hash = it->first;
                    reaction_id = it->second[i].second;

                    std::size_t pos = hash.find(".");
                    site_one = std::optional<int>(stoi(hash.substr(0, pos)));
                    site_two = std::optional<int>(stoi(hash.substr(pos+1)));
                    if(site_one < site_two) {
                        assert(false);
                    }
                    isFound = true;
                    break;
                }
            }
            it++;
        }
    }
    
    dt = - std::log(r2) / propensity_sum;


    if(isFound) {
        return std::optional<LatticeEvent> ( LatticeEvent {.index = reaction_id, .dt = dt,
                                                           .site_one = site_one, .site_two = site_two});
    }
    else {

        return std::optional<LatticeEvent> ();
    }
        
}

/* ---------------------------------------------------------------------- */

std::string LatSolver::make_string(int site_one, int site_two) {

    return (site_one > site_two) ? std::to_string(site_one) + "." + 
    std::to_string(site_two) : std::to_string(site_two) + "." + std::to_string(site_one);

} // make_string


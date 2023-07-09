#ifndef RNMC_LATTICE_SOLVER_H
#define RNMC_LATTICE_SOLVER_H

#include "../core/sampler.h"
#include "../core/RNMC_types.h"

#include <vector>
#include <unordered_map>
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

class LatticeSolver {
public:
    LatticeSolver() : sampler(Sampler(0)) {}; // defualt constructor
    LatticeSolver(unsigned long int seed, std::vector<double> &&initial_propensities);
    LatticeSolver(unsigned long int seed, std::vector<double> &initial_propensities);

    void update(Update update);
    void update(std::vector<Update> updates);

    void update(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                std::vector< std::pair<double, int> > > &props);
    void update(std::vector<LatticeUpdate> lattice_updates, 
                std::unordered_map<std::string,                     
                std::vector< std::pair<double, int> > > &props);

    std::optional<LatticeEvent> event_lattice(std::unordered_map<std::string,                     
                                            std::vector< std::pair<double, int> > > &props);

    std::string make_string(int site_one, int site_two);

    long double propensity_sum;
    int number_of_active_indices;               // end simulation of no sites with non zero propensity   

    std::vector<double> propensities;                   // Gillepsie propensities 

private:
    Sampler sampler;           
};

#endif
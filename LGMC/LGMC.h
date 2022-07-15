
#ifndef LGMC_H
#define LGMC_H

#include "stdio.h"
#include "memory.h"
#include "lattice.h"
#include "../GMC/reaction_network.h"
#include "../core/sampler.h"


namespace LGMC_NS {

class LGMC {
    public: 
        LGMC(int argc, char **argv);
        
        ~LGMC();
        
        void print_usage();

        void run();


        /* -------------------------------- Lattice Updates ----------------------------- */

        bool update_state(int next_reaction, int site_one, int site_two);

        void clear_site(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site, std::optional<int> ignore_neighbor);

        void clear_site_helper(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site_one, int site_two);

        void relevant_react(std::function<void(LatticeUpdate lattice_update)>
                            update_function, int site, std::optional<int> ignore_neighbor);

        void compute_propensitycompute_propensity(int site_one, int site_two, int num_one, 
                                                  int num_two, int react_id)

        void update_propensities(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        std::function<void(LatticeUpdate lattice_update)> update_function,
                        int next_reaction);

        void update_electrochemical();

        std::string make_string(int site_one, int site_two);

        /* -------------------------------- Gillespie Updates ----------------------------- */
        void update_state(std::vector<int> &state,
        int reaction_index);

        void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction
        );


    private:                                                          

        int prop_sum;                               // running total of propensities
                                    
        LatticeReactionNetwork *react_net;                   // pointer to gillespie reaction network

        Memory *memory;                            // memory allocation functions
        Lattice *lattice;                          // lattice for SEI
        Sampler sampler;

};

}

#endif


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
        
        void update_propensity(int site_one, 
        int site_two, LatticeReaction *reaction, int react_id, int num_one, int num_two);

        void update(std::optional<int> site_one, std::optional<int> site_two, int reaction_id);

        void relevant_react(int site);

        void clear_site(int site, std::optional<int> ignore_neighbor);

        void clear_site_helper(int site_one, int site_two);

        void run();

        void update_one_lattice_site(int site, LatticeReaction *reaction);

        void update_two_lattice_sites(int site_one, int site_two, LatticeReaction *reaction);

        int sum_row(std::string hash);

        std::string make_string(int site_one, int site_two);

        void update_site(int site, std::optional<int> ignore_neighbor);

        void update_electrochemical();


    public: 
    
        std::unordered_map<std::string,             // lattice propensities as site neighbor pair
            std::vector< std::pair<double, int> > > props;          
                                                      

        int prop_sum;                               // running total of propensities
                                    
        LatticeReactionNetwork *react_net;                   // pointer to gillespie reaction network

        Memory *memory;                            // memory allocation functions
        Lattice *lattice;                          // lattice for SEI
        Sampler sampler;

};

}

#endif

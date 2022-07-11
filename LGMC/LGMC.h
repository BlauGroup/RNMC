
#ifndef LGMC_H
#define LGMC_H

//#include "mpi.h"
#include "stdio.h"
#include "memory.h"
#include "lattice.h"
#include "../core/dispatcher.h"
#include "../GMC/reaction_network.h"
#include "../GMC/sql_types.h"


namespace LGMC_NS {

class LGMC {
    public: 
        Memory *memory;                            // memory allocation functions
        Lattice *lattice;                          // lattice for SEI


        LGMC(int argc, char **argv);
        
        ~LGMC();
        
        void print_usage();
        
        void update_propensity(int site_one, 
        int site_two, Reaction &reaction, int react_id);

        void update(int site_one, int site_two);

        void relevant_react(int site);

        void clear_site(int site);

        void run();

        int sum_row(std::string hash);

        std::string make_string(int site_one, int site_two);


    private: 
    
        std::unordered_map<std::string,             // lattice propensities as site neighbor pair
            std::vector< std::pair<double, int> > > props;          
                                                    

        std::vector<std::vector<int>>               // dependency list without products  
        lat_dependents;         

        int prop_sum;                               // running total of propensities
                                    
        ReactionNetwork *react_net;                   // pointer to gillespie reaction network

};

}

#endif


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
        LGMC();
        
        ~LGMC();
        
        void print_usage();

        /* -------------------------------- Lattice Updates ----------------------------- */

        bool update_state(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int next_reaction, int site_one, int site_two, double &prop_sum);

        void clear_site(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site, std::optional<int> ignore_neighbor, double &prop_sum);

        void clear_site_helper(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site_one, int site_two, double &prop_sum);

        void relevant_react(std::function<void(LatticeUpdate lattice_update)> update_function, 
                            std::vector<int> &state, int site, std::optional<int> ignore_neighbor);

        double compute_propensity(int site_one, int site_two, int num_one, 
                                int num_two, int react_id);

        void update_propensities(std::vector<int> &state,
                        std::function<void(LatticeUpdate lattice_update)> 
                        update_function, int next_reaction, int site_one, int site_two);

        std::string make_string(int site_one, int site_two);

        double sum_row(std::string hash, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props);

        /* for EVAN */
        void update_electrochemical();

        void update_adsorption();

        /* -------------------------------- Gillespie Updates ----------------------------- */
        void update_state(std::vector<int> &state, int reaction_index);

        void update_propensities(std::function<void(Update update)> update_function,
                                 std::vector<int> &state, int next_reaction);

        // convert a history element as found a simulation to history
        // to a SQL type.
        TrajectoriesSql history_element_to_sql(
            int seed,
            HistoryElement history_element);

    private:                                                          
                           
        LatticeReactionNetwork *react_net;                   // pointer to gillespie reaction network

        Memory *memory;                            // memory allocation functions
        Lattice *lattice;                          // lattice for SEI
        Sampler sampler;

};

}

#endif

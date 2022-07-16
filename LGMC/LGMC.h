
#ifndef LGMC_H
#define LGMC_H

#include "stdio.h"
#include "memory.h"
#include "lattice.h"
#include "../GMC/reaction_network.h"
#include "../core/sampler.h"
#include "../core/simulation.h"
#include "LatSolver.h"

struct LGMCParameters {
    float latconst;                               
    float boxxlo,boxxhi,boxylo,                   
          boxyhi,boxzlo,boxzhi;                       
    bool xperiodic, yperiodic, zperiodic;
    float temperature;
    float potential;
};

namespace LGMC_NS {

class LGMC {
    public: 
        LGMC(SqlConnection &reaction_network_database, SqlConnection &initial_state_database, 
            LGMCParameters parameters);
        
        ~LGMC();

        /* -------------------------------- Updates LGMC ----------------------------- */

        void update_state(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        std::vector<int> &state, int next_reaction, 
                        std::optional<int> site_one, std::optional<int> site_two, double prop_sum);

        void update_propensities(std::vector<int> &state, this->update_function, 
                                        lattice_update_function, int next_reaction, 
                                        std::optional<int> site_one, std::optional<int> site_two);

        void update_adsorption(); // TODO

        /* -------------------------------- Updates Lattice ----------------------------- */

        bool update_state(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int next_reaction, int site_one, int site_two, double prop_sum);

        void clear_site(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site, std::optional<int> ignore_neighbor, double prop_sum);

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

        /* -------------------------------------------------------------------------------- */
        // convert a history element as found a simulation to history
        // to a SQL type.
        TrajectoriesSql history_element_to_sql(
            int seed,
            HistoryElement history_element);
        
        std::vector<double> initial_propensities;
        std::vector<int> initial_state;

    private:                                                          
                           
        LatticeReactionNetwork *react_net;                   // pointer to gillespie reaction network

        Memory *memory;                            // memory allocation functions
        Lattice *lattice;                          // lattice for SEI
        Sampler sampler;

};

}

#endif

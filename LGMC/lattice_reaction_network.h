
#ifndef LGMC_H
#define LGMC_H

#include "stdio.h"
#include "memory.h"
#include "lattice.h"
#include "../core/sql.h"
#include "../core/sampler.h"
#include "../core/simulation.h"
#include "../core/solvers.h"
#include "LatSolver.h"
#include "../GMC/sql_types.h"


enum Phase {LATTICE, SOLUTION};
enum Type {ADSORPTION, DESORPTION, HOMOGENEOUS_ELYTE, HOMOGENEOUS_SOLID, DIFFUSION, OXIDATION, REDUCTION};
enum ChargeTransferStyle {MARCUS, BUTLER_VOLMER};

struct LGMCParameters {
    float latconst;                               
    float boxxlo,boxxhi,boxylo,                   
          boxyhi,boxzlo,boxzhi;                       
    bool xperiodic, yperiodic, zperiodic;
    float temperature;
    float g_e;
    bool is_add_sites;
    ChargeTransferStyle charge_transfer_style;
};


struct LatticeReaction {
    // we assume that each reaction has zero, one or two reactants
    uint8_t number_of_reactants;
    uint8_t number_of_products;

    int reactants[2];
    int products[2];

    Phase phase_reactants[2];
    Phase phase_products[2];

    double dG;
    double reorganization;
    double rate;

    Type type;
};

namespace LGMC_NS {

class LatticeReactionNetwork {
    public: 
        LatticeReactionNetwork(SqlConnection &reaction_network_database, SqlConnection &initial_state_database, 
            LGMCParameters parameters);
        
        ~LatticeReactionNetwork();

        /* -------------------------------- Updates Global ----------------------------- */

        void update_state(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        std::vector<int> &state, int next_reaction, 
                        std::optional<int> site_one, std::optional<int> site_two, 
                        double &prop_sum, int &active_indices);

        void update_propensities(Lattice *lattice, std::vector<int> &state, std::function<void(Update update)> update_function, 
                                        std::function<void(LatticeUpdate lattice_update)> lattice_update_function, 
                                        int next_reaction, std::optional<int> site_one, std::optional<int> site_two);

        void update_adsorp_state(Lattice *lattice, std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props,
                                double &prop_sum, int &active_indices); 

        void update_adsorp_props(Lattice *lattice, std::function<void(LatticeUpdate lattice_update)> lattice_update_function, std::vector<int> &state);

        /* -------------------------------- Updates Lattice ----------------------------- */

        bool update_state(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int next_reaction, int site_one, int site_two, 
                        double &prop_sum, int &active_indices);

        void clear_site(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site, std::optional<int> ignore_neighbor, 
                        double &prop_sum, int &active_indices);

        void clear_site_helper(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site_one, int site_two, double &prop_sum, 
                        int &active_indices);

        void relevant_react(Lattice *lattice, std::function<void(LatticeUpdate lattice_update)> update_function,
                                            int site, std::optional<int> ignore_neighbor);

        double compute_propensity(int num_one, int num_two, int react_id, Lattice *lattice, int site_id = 0);

        bool update_propensities(Lattice *lattice, 
                        std::function<void(LatticeUpdate lattice_update)> 
                        update_function, int next_reaction, int site_one, int site_two);

        std::string make_string(int site_one, int site_two);

        double sum_row(std::string hash, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props);
        
        /* -------------------------- Updates Reaction Network ----------------------------- */

        void init_reaction_network(SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database, Lattice *lattice);

        void fill_reactions(SqlConnection &reaction_network_database);

        void compute_dependents();

        double compute_propensity(std::vector<int> &state, int reaction_index, Lattice *lattice);

        void update_propensities(std::function<void(Update update)> update_function,
                                std::vector<int> &state, int next_reaction, Lattice *lattice);
        
        // for compatibility 
        void update_propensities(std::function<void(Update update)> update_function,
                                std::vector<int> &state, int next_reaction) {assert(false);};

        void update_state(std::vector<int> &state, int reaction_index);

        /* -------------------------------------------------------------------------------- */
        
        std::vector<double> initial_propensities;
        std::vector<int> initial_state;
        Lattice *initial_lattice;                          // lattice for SEI

        // convert a history element as found a simulation to history
        // to a SQL type.
        TrajectoriesSql history_element_to_sql(
            int seed,
            HistoryElement history_element);

    private:                                                          

        Sampler sampler;
        std::vector<LatticeReaction> reactions;

        double factor_two; // rate modifier for reactions with two reactants
        double factor_duplicate; // rate modifier for reactions of form A + A -> ...
        std::vector<std::vector<int>> dependents;
        
        bool is_add_sites;
        float temperature;
        float g_e;

};

}

#endif

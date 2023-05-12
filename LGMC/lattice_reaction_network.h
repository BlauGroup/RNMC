#ifndef RNMC_LATTICE_REACTION_NETWORK_H
#define RNMC_LATTICE_REACTION_NETWORK_H

#include "memory.h"
#include "lattice.h"
#include "lattice_solver.h"
#include "sql_types.h"

#include "../core/sql.h"
#include "../core/sampler.h"
#include "../core/sql_types.h"
#include "../core/queues.h"


#include <list>
#include <vector>
#include <string>
#include <assert.h>

const int SITE_SELF_REACTION = -3;
const int SITE_HOMOGENEOUS = -2;

const double KB = 8.6173e-5;           // In eV/K
const double PLANCK = 4.1357e-15;      // In eV s


enum Phase {LATTICE, SOLUTION, NONE};
enum Type {ADSORPTION, DESORPTION, HOMOGENEOUS_ELYTE, HOMOGENEOUS_SOLID, DIFFUSION, OXIDATION, REDUCTION};
enum ChargeTransferStyle {MARCUS, BUTLER_VOLMER};

struct LatticeState {
    std::vector<int> homogeneous;
    Lattice *lattice;
};

struct LatticeParameters {
    float latconst;                               
    float boxxlo,boxxhi,boxylo,                   
          boxyhi,boxzlo,boxzhi;                       
    bool xperiodic, yperiodic, zperiodic;
    float temperature;
    float g_e;
    bool is_add_sites;
    ChargeTransferStyle charge_transfer_style;
    std::string lattice_fill;
};


struct LatticeReaction {
    // we assume that each reaction has zero, one or two reactants
    uint8_t number_of_reactants;
    uint8_t number_of_products;

    int reactants[2];
    int products[2];

    Phase phase_reactants[2];
    Phase phase_products[2];

    // Basic thermodynamics/kinetics
    double dG; // eV
    double prefactor; // Hz (s^-1)
    double rate; // Units depend on reaction molecularity

    // Charge transfer
    double electron_tunneling_coefficient; // A^-1
    
    // Marcus theory
    double reorganization_energy; // eV

    // Butler-Volmer theory
    double charge_transfer_coefficient; // Unitless

    Type type;
};

class LatticeReactionNetwork {
    public: 
        LatticeReactionNetwork(SqlConnection &reaction_network_database, SqlConnection &initial_state_database, 
            LatticeParameters parameters);
        
        ~LatticeReactionNetwork();

        /* -------------------------------- Updates Global ----------------------------- */

        void update_state(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        std::vector<int> &state, int next_reaction, 
                        std::optional<int> site_one, std::optional<int> site_two, 
                        long double &prop_sum, int &active_indices);

        void update_propensities(Lattice *lattice, std::vector<int> &state, std::function<void(Update update)> update_function, 
                                        std::function<void(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props)> lattice_update_function, 
                                        int next_reaction, std::optional<int> site_one, std::optional<int> site_two, 
                                        std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props);

        void update_adsorp_state(Lattice *lattice, std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props,
                                long double &prop_sum, int &active_indices); 

        void update_adsorp_props(Lattice *lattice, std::function<void(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props)> lattice_update_function, std::vector<int> &state, 
                        std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props);

        /* -------------------------------- Updates Lattice ----------------------------- */

        bool update_state(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int next_reaction, int site_one, int site_two, 
                        long double &prop_sum, int &active_indices);

        void clear_site(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site, std::optional<int> ignore_neighbor, 
                        long double &prop_sum, int &active_indices);

        void clear_site_helper(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site_one, int site_two, long double &prop_sum, 
                        int &active_indices);

        void relevant_react(Lattice *lattice, std::function<void(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props)> update_function,
                                            int site, std::optional<int> ignore_neighbor,
                                            std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props);

        double compute_propensity(int num_one, int num_two, int react_id, Lattice *lattice, int site_id = 0);

        bool update_propensities(Lattice *lattice, 
                        std::function<void(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props)> 
                        update_function, int next_reaction, int site_one, int site_two, 
                        std::unordered_map<std::string,std::vector< std::pair<double, int> > > &props);

        std::string make_string(int site_one, int site_two);
        std::string make_string(std::vector<int> vec);

        double sum_row(std::string hash, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props);
        
        void update_all_propensities(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        long double &prop_sum, int &active_indices,
                        std::function<void(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props)> 
                        update_function);

        /* -------------------------- Updates Reaction Network ----------------------------- */

        void init_reaction_network(SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database, Lattice *lattice);

        void fill_reactions(SqlConnection &reaction_network_database);

        void compute_dependents();

        double compute_propensity(std::vector<int> &state, int reaction_index, Lattice *lattice);

        void update_propensities(std::function<void(Update update)> update_function,
                                std::vector<int> &state, int next_reaction, Lattice *lattice);
        
        void update_state(std::vector<int> &state, int reaction_index);

        void compute_initial_propensities(std::vector<int> state, Lattice *lattice);

        /* -------------------------------------------------------------------------------- */
        
        std::vector<double> initial_propensities;
        std::vector<int> initial_state;
        Lattice *initial_lattice;                          // lattice for SEI

        // convert a history element as found a simulation to history
        // to a SQL type.
        LatticeWriteTrajectoriesSql history_element_to_sql(
            int seed,
            LatticeTrajectoryHistoryElement history_element);

        LatticeWriteStateSql state_history_element_to_sql
            (int seed, LatticeStateHistoryElement history_element);

        LatticeWriteCutoffSql cutoff_history_element_to_sql(
            int seed,
            LatticeCutoffHistoryElement cutoff_history_element);

        void checkpoint(SqlReader<LatticeReadStateSql> state_reader, 
                                        SqlReader<LatticeReadCutoffSql> cutoff_reader, 
                                        SqlReader<LatticeReadTrajectoriesSql> trajectory_reader, 
                                        std::map<int, LatticeState> &temp_seed_state_map, 
                                        std::map<int, int> &temp_seed_step_map, 
                                        SeedQueue &temp_seed_queue, 
                                        std::map<int, double> &temp_seed_time_map, 
                                        LatticeReactionNetwork &model);

        void store_checkpoint(std::vector<LatticeStateHistoryElement> &state_packet,
        LatticeState &state, LatticeReactionNetwork &lattice_reaction_network, unsigned long int &seed, 
        int step, double time, std::vector<LatticeCutoffHistoryElement> &cutoff_packet);


        double get_butler_volmer_rate_coefficient(double base_dg, double prefactor, double charge_transfer_coefficient,
                                          double electron_tunneling_coefficient, double e_free, double distance,
                                          double temperature, bool reduction);
        

        double get_marcus_rate_coefficient(double base_dg, double prefactor, double reorganization_energy,
                                   double electron_tunneling_coefficient, double e_free, double distance,
                                   double temperature, bool reduction);

        void print_state_propensities(long double propensity_sum,
                        std::vector<double> &propensities,
                        std::vector<int> &state, std::string filename);

        int szudzik(int a, int b);

        int combine(int i, int j, int k);

        std::unordered_map<int, std::tuple<uint32_t,uint32_t,uint32_t>>
        szudzik_mapping(int i_max, int j_max, int k_max);
    
    private:                                                          

        Sampler sampler;
        std::vector<LatticeReaction> reactions;
        std::unordered_map<int,int> species_size;           // key: species ID, value: size of the species

        double factor_two; // rate modifier for reactions with two reactants
        double factor_duplicate; // rate modifier for reactions of form A + A -> ...
        std::vector<std::vector<int>> dependents;
        
        ChargeTransferStyle charge_transfer_style;

        bool is_add_sites;
        float temperature;
        float g_e;

};

#endif
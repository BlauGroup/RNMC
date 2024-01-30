#ifndef RNMC_LATTICE_REACTION_NETWORK_H
#define RNMC_LATTICE_REACTION_NETWORK_H

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
#include <fstream>
#include <functional>
#include <memory>

const int SITE_SELF_REACTION = -3;
const int SITE_HOMOGENEOUS = -2;

const double KB = 8.6173e-5;           // In eV/K
const double PLANCK = 4.1357e-15;      // In eV s

enum Phase {LATTICE, SOLUTION, NONE};
enum Type {ADSORPTION, DESORPTION, HOMOGENEOUS_SOLUTION, 
           HOMOGENEOUS_LATTICE, DIFFUSION, OXIDATION, REDUCTION};
enum ChargeTransferStyle {MARCUS, BUTLER_VOLMER};

class LatticeState {

public:
    std::vector<int> homogeneous;
    std::unique_ptr<Lattice> lattice;
    
    LatticeState(); // default constructor
    LatticeState(std::vector<int> homogeneous_in, std::unique_ptr<Lattice> lattice_in);
    LatticeState(const LatticeState & lattice_in);
};

struct LatticeParameters {
    float latconst;                               
    float boxxhi, boxyhi, boxzhi;                       
    float temperature;
    float g_e;
    bool is_add_sites;
    ChargeTransferStyle charge_transfer_style;
    bool isCheckpoint;
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
    LatticeReactionNetwork(); // default constructor
    
    LatticeReactionNetwork(SqlConnection &reaction_network_database, 
                        SqlConnection &initial_state_database, 
                        LatticeParameters parameters);
    
   // ~LatticeReactionNetwork();

    /* -------------------------------- Updates Global ----------------------------- */

    void update_state(std::unique_ptr<Lattice> &lattice, std::unordered_map<std::string,                     
                    std::vector< std::pair<double, int> > > &props,
                    std::vector<int> &state, int next_reaction, 
                    std::optional<int> site_one, std::optional<int> site_two, 
                    long double &prop_sum, int &active_indices, bool &flip_sites);

    void update_propensities(std::unique_ptr<Lattice> &lattice, std::vector<int> &state, 
                            std::function<void(Update update)> update_function, 
                            std::function<void(LatticeUpdate lattice_update, std::unordered_map<std::string,
                            std::vector< std::pair<double, int> > > &props)> lattice_update_function, 
                            int next_reaction, std::optional<int> site_one, std::optional<int> site_two, 
                            std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props);

    void update_adsorp_state(std::unique_ptr<Lattice> &lattice, std::unordered_map<std::string, 
                            std::vector< std::pair<double, int> > > &props,
                            long double &prop_sum, int &active_indices); 

    void update_adsorp_props(std::unique_ptr<Lattice> &lattice, std::function<void(LatticeUpdate lattice_update, 
                            std::unordered_map<std::string,                     
                            std::vector< std::pair<double, int> > > &props)> lattice_update_function, 
                            std::vector<int> &state, 
                            std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props);

    /* -------------------------------- Updates Lattice ----------------------------- */

    bool update_state_lattice(std::unique_ptr<Lattice> &lattice, std::unordered_map<std::string,                     
                            std::vector< std::pair<double, int> > > &props, 
                            int next_reaction, int site_one, int site_two, 
                            long double &prop_sum, int &active_indices, bool &flip_sites);

    void clear_site(std::unique_ptr<Lattice> &lattice, std::unordered_map<std::string,                     
                    std::vector< std::pair<double, int> > > &props,
                    int site, std::optional<int> ignore_neighbor, 
                    long double &prop_sum, int &active_indices);

    void clear_site_helper(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site_one, int site_two, long double &prop_sum, 
                        int &active_indices);

    void relevant_react(std::unique_ptr<Lattice> &lattice, std::function<void(LatticeUpdate lattice_update, 
                    std::unordered_map<std::string,                     
                    std::vector< std::pair<double, int> > > &props)> update_function,
                    int site, std::optional<int> ignore_neighbor,
                    std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props);

    double compute_propensity(int num_one, int num_two, int react_id, std::unique_ptr<Lattice> &lattice, int site_id = 0);

    bool update_propensities(std::unique_ptr<Lattice> &lattice, 
                        std::function<void(LatticeUpdate lattice_update, 
                        std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props)> 
                        update_function, int next_reaction, int site_one, int site_two, 
                        std::unordered_map<std::string,std::vector< std::pair<double, int> > > &props);

    std::string make_string(int site_one, int site_two);
    std::string make_string(std::vector<int> vec);

    double sum_row(std::string hash, std::unordered_map<std::string,                     
                    std::vector< std::pair<double, int> > > &props);
    
    void update_all_propensities(std::unique_ptr<Lattice> &lattice, std::unordered_map<std::string,                     
                            std::vector< std::pair<double, int> > > &props, 
                            long double &prop_sum, int &active_indices,
                            std::function<void(LatticeUpdate lattice_update, 
                            std::unordered_map<std::string,                     
                            std::vector< std::pair<double, int> > > &props)> 
                            update_function);

    /* -------------------------- Updates Reaction Network ----------------------------- */

    void init_reaction_network(SqlConnection &reaction_network_database,
    SqlConnection &initial_state_database,
    LatticeParameters parameters);

    void fill_reactions(SqlConnection &reaction_network_database);

    void compute_dependents();

    double compute_propensity(std::vector<int> &state, int reaction_index, std::unique_ptr<Lattice> &lattice );

    void update_propensities(std::function<void(Update update)> update_function,
                            std::vector<int> &state, int next_reaction, std::unique_ptr<Lattice> &lattice);
    
    void update_state_solution(std::vector<int> &state, int reaction_index);

    void compute_initial_propensities(std::vector<int> state, std::unique_ptr<Lattice> &lattice);

    /* -------------------------------------------------------------------------------- */                       

    // convert a history element as found a simulation to history
    // to a SQL type.
    LatticeWriteTrajectoriesSql history_element_to_sql(int seed,
        LatticeTrajectoryHistoryElement history_element);

    LatticeWriteStateSql state_history_element_to_sql(int seed, 
        LatticeStateHistoryElement history_element);

    LatticeWriteCutoffSql cutoff_history_element_to_sql(int seed,
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
                        LatticeState &state, unsigned long int &seed, 
                        int step, double time, 
                        std::vector<LatticeCutoffHistoryElement> &cutoff_packet);


    double get_butler_volmer_rate_coefficient(double base_dg, double prefactor, 
                                            double charge_transfer_coefficient,
                                            double electron_tunneling_coefficient, 
                                            double e_free, double distance,
                                            double temperature, bool reduction);
    
    double get_marcus_rate_coefficient(double base_dg, double prefactor, 
                                    double reorganization_energy,
                                    double electron_tunneling_coefficient, 
                                    double e_free, double distance,
                                    double temperature, bool reduction);

    void print_state_propensities(long double propensity_sum,
                                std::vector<double> &propensities,
                                std::vector<int> &state, std::string filename);

    int szudzik(int a, int b);

    int combine(int i, int j, int k);

    std::unordered_map<int, std::tuple<uint32_t,uint32_t,uint32_t>>
        szudzik_mapping(int i_max, int j_max, int k_max);

    std::vector<double> initial_propensities;
    LatticeState initial_state; 

    std::vector<LatticeReaction> reactions;
    std::vector<std::vector<int>> dependents;
    bool isCheckpoint;
    
private:                                                          

    Sampler sampler;
    std::unordered_map<int,int> species_size;           // key: species ID, value: size of the species

    double factor_two; // rate modifier for reactions with two reactants
    double factor_duplicate; // rate modifier for reactions of form A + A -> ...
    
    ChargeTransferStyle charge_transfer_style;

    bool is_add_sites;
    float temperature;
    float g_e;
};

#endif
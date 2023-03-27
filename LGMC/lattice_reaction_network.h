#pragma once
#include "stdio.h"
#include "memory.h"
#include "lattice.h"
#include "../core/sql.h"
#include "../core/sampler.h"
#include "../core/sql_types.h"
#include "../core/solvers.h"
#include "LatSolver.h"
#include "../LGMC/sql_types.h"

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
        
        /* -------------------------- Updates Reaction Network ----------------------------- */

        void init_reaction_network(SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database, Lattice *lattice);

        void fill_reactions(SqlConnection &reaction_network_database);

        void compute_dependents();

        double compute_propensity(std::vector<int> &state, int reaction_index, Lattice *lattice);

        void update_propensities(std::function<void(Update update)> update_function,
                                std::vector<int> &state, int next_reaction, Lattice *lattice);
        
        void update_state(std::vector<int> &state, int reaction_index);

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

        WriteCutoffSql cutoff_history_element_to_sql(
            int seed,
            CutoffHistoryElement cutoff_history_element);

        bool read_state(SqlReader<LatticeReadStateSql> state_reader, 
                    std::map<int, LatticeState> &temp_seed_state_map,
                LatticeReactionNetwork &lattice_reaction_network,
                SeedQueue &temp_seed_queue);

        void read_trajectories(SqlReader<LatticeReadTrajectoriesSql> trajectory_reader, 
                           std::map<int, LatticeState> &temp_seed_state_map, 
                           std::map<int, int> &temp_seed_step_map, 
                           std::map<int, double> &temp_seed_time_map,
                           LatticeReactionNetwork &lattice_reaction_network);

        void store_state_history(std::vector<LatticeStateHistoryElement> &state_packet,
        LatticeState &state, LatticeReactionNetwork &lattice_reaction_network, unsigned long int &seed);
           

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


double get_butler_volmer_rate_coefficient(double base_dg, double prefactor, double charge_transfer_coefficient,
                                          double electron_tunneling_coefficient, double e_free, double distance,
                                          double temperature, bool reduction) {

    double dg, kappa;

    if (reduction) {
        dg = base_dg - e_free;
    }
    else {
        dg = base_dg + e_free;
    }

    kappa = std::exp(-1 * electron_tunneling_coefficient * distance);
    return kappa * prefactor * std::exp(-1 * charge_transfer_coefficient * dg / (KB * temperature));

}


double get_marcus_rate_coefficient(double base_dg, double prefactor, double reorganization_energy,
                                   double electron_tunneling_coefficient, double e_free, double distance,
                                   double temperature, bool reduction) {

    double dg, dg_barrier, squared, kappa;

    if (reduction) {
        dg = base_dg - e_free;
    }
    else {
        dg = base_dg + e_free;
    }

    squared = 1 + dg / reorganization_energy;
    dg_barrier = reorganization_energy / 4 * squared * squared;
    kappa = std::exp(-1 * electron_tunneling_coefficient * distance);

    if (dg_barrier < 0) {
        return kappa * prefactor;
    } else {
        return kappa * prefactor * std::exp(-1 * dg_barrier / (KB * temperature));
    }
}


LatticeReactionNetwork::LatticeReactionNetwork(SqlConnection &reaction_network_database, SqlConnection &initial_state_database, 
            LatticeParameters parameters) : sampler (Sampler(0)) {

    // create lattice
    initial_lattice = new Lattice(parameters.latconst, parameters.boxxlo, parameters.boxxhi, parameters.boxylo,
                    parameters.boxyhi, parameters.boxzlo, parameters.boxzhi, parameters.xperiodic, parameters.yperiodic, parameters.zperiodic);


    if(parameters.lattice_fill.compare("none") != 0) {
        // fill lattice
        initial_lattice->fill(parameters.lattice_fill);
    }
    init_reaction_network(reaction_network_database, initial_state_database, initial_lattice);

    is_add_sites = parameters.is_add_sites;
    temperature = parameters.temperature;
    g_e = parameters.g_e;
    charge_transfer_style = parameters.charge_transfer_style;

} // LatticeReactionNetwork()

/* ---------------------------------------------------------------------- */

LatticeReactionNetwork::~LatticeReactionNetwork()
{
    delete initial_lattice;
} // ~LatticeReactionNetwork()

/* ---------------------------------------------------------------------- 
    Global Updates
---------------------------------------------------------------------- */

void LatticeReactionNetwork::update_state(Lattice *lattice,
                        std::unordered_map< std::string, std::vector< std::pair< double, int > > > &props,
                        std::vector<int> &state, int next_reaction, 
                        std::optional<int> site_one, std::optional<int> site_two, 
                        long double &prop_sum, int &active_indices) {
    
    if(site_one) {
        // update lattice state
        bool update_gillepsie = update_state(lattice, props, next_reaction, site_one.value(), 
                                            site_two.value(), prop_sum, active_indices);
        
        if(update_gillepsie) {
            update_state(state, next_reaction);
        }
    }
    
    else {
        // homogeneous event happens
        update_state(state, next_reaction);

    }

    update_adsorp_state(lattice, props, prop_sum, active_indices);

} // update_state

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::update_propensities(Lattice *lattice, std::vector<int> &state, std::function<void(Update update)> update_function, 
                                        std::function<void(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props)> lattice_update_function, 
                                        int next_reaction, std::optional<int> site_one, std::optional<int> site_two, 
                                        std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props) {

    if(site_one) {
        // update lattice state
        bool update_gillepsie = update_propensities(lattice, lattice_update_function, next_reaction, site_one.value(), site_two.value(), props);
        
        if(update_gillepsie) {
            update_propensities(update_function, state, next_reaction, lattice);
        }
    }
    
    else {
        // homoegenous event happens
        update_propensities(update_function, state, next_reaction, lattice);

    }

    update_adsorp_props(lattice, lattice_update_function, state, props);

}

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::update_adsorp_state(Lattice *lattice, std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props,
                                long double &prop_sum, int &active_indices) {
    // update only sites on the edge 
    for(auto it = lattice->edges.begin(); it != lattice->edges.end(); it++) {
        int site = it->first;
        if(lattice->edges[site] == 'a') {
             assert(lattice->sites[site].species == SPECIES_EMPTY);
             clear_site_helper(props, site, SITE_HOMOGENEOUS, prop_sum, active_indices);
        }
       
    }

} // update_adsorp_state

void LatticeReactionNetwork::update_adsorp_props(Lattice *lattice, std::function<void(LatticeUpdate lattice_update, 
std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props)> lattice_update_function,
                                                 std::vector<int> &state, std::unordered_map<std::string,                     
                std::vector< std::pair<double, int> > > &props) {

    // update only sites on the edge 
    for(auto it = lattice->edges.begin(); it != lattice->edges.end(); it++) {
        int site = it->first;

        if(lattice->edges[site] == 'a') {
            // find relevant adsorption reactions
            // the species on adsoprtion sites will always be empty so dependence
            // graph will point to all adsorption reactions
            std::vector<int> &potential_reactions = dependents[lattice->sites[site].species]; 

            for(unsigned long int i = 0; i < static_cast<int> (potential_reactions.size()); i++ ) {

                unsigned long int reaction_id = potential_reactions[i];
                LatticeReaction reaction = reactions[reaction_id]; 

                if(reaction.type == Type::ADSORPTION) {

                    int other_reactant_id = (reaction.reactants[0] == lattice->sites[site].species) ? 1 : 0;

                    if(state[reaction.reactants[other_reactant_id]] != 0) {
                        
                        double new_propensity = compute_propensity(1, state[reaction.reactants[other_reactant_id]], reaction_id, lattice);

                        lattice_update_function(LatticeUpdate {
                            .index = reaction_id,
                            .propensity = new_propensity,
                            .site_one = site,
                            .site_two = SITE_HOMOGENEOUS}, props); 

                    }
                    
                    
                }

            }

        }

    }

}
/* ---------------------------------------------------------------------- 
---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- 
    Only calls this function if necessary reactants are on lattice sites
---------------------------------------------------------------------- */

double LatticeReactionNetwork::compute_propensity(int num_one, int num_two, 
                                int react_id, Lattice *lattice, int site_id) {   

    LatticeReaction reaction = reactions[react_id]; 

    double p, k;

    if(reaction.type == Type::REDUCTION || reaction.type == Type::OXIDATION) {
        bool reduction_in = (reaction.type == Type::REDUCTION) ? true : false;

        assert(charge_transfer_style == ChargeTransferStyle::MARCUS || charge_transfer_style == ChargeTransferStyle::BUTLER_VOLMER);

        // TODO: make general for all types of periodicity (distance = .z)
        if (charge_transfer_style == ChargeTransferStyle::MARCUS) {
            k = get_marcus_rate_coefficient(reaction.dG, reaction.prefactor, reaction.reorganization_energy,
                                            reaction.electron_tunneling_coefficient, g_e, lattice->sites[site_id].z, temperature,
                                            reduction_in);
        } else {
            k = get_butler_volmer_rate_coefficient(reaction.dG, reaction.prefactor, reaction.charge_transfer_coefficient,
                                                   reaction.electron_tunneling_coefficient, g_e, lattice->sites[site_id].z, temperature,
                                                   reduction_in);    
        }
    } else {
        k = reaction.rate;
    }

    // one reactant
    if (reaction.number_of_reactants == 1)
        p = num_one * k;


    // two reactants
    else {
        if (reaction.reactants[0] == reaction.reactants[1])
            p = factor_duplicate
                * factor_two
                * num_one
                * (num_one - 1)
                * k;

        else
            p = factor_two
                * num_one
                * num_two
                * k;
    }

    assert(p > 0);

    return p;
    
} // compute_propensity() 

/* ---------------------------------------------------------------------- */

bool LatticeReactionNetwork::update_state(Lattice *lattice, 
                        std::unordered_map<std::string, std::vector< std::pair< double, int > > > &props, 
                        int next_reaction, int site_one, int site_two, 
                        long double &prop_sum, int &active_indices) {

    LatticeReaction reaction = reactions[next_reaction]; 

    if(reaction.type == Type::ADSORPTION) {

        assert(lattice->edges.find(site_one) != lattice->edges.end());
        assert(lattice->edges[site_one] == 'a');
        assert(lattice->sites[site_one].species == SPECIES_EMPTY);
        assert(site_two == SITE_HOMOGENEOUS);
        // TODO: The below might be a slow operation

        // update site
        lattice->sites[site_one].species = reaction.products[0];
        clear_site(lattice, props, site_one, std::optional<int> (), prop_sum, active_indices);
        clear_site_helper(props, site_one, SITE_HOMOGENEOUS, prop_sum, active_indices);

        if(is_add_sites) {
            // TODO: make general for all types of periodicity (add_site coordinates)
            lattice->add_site(lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1,
            lattice->sites[site_one].x, lattice->sites[site_one].y, lattice->sites[site_one].z + lattice->get_latconst(), 
            true, true, true);
            
            std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1};
            int site_new = lattice->loc_map[key];

            clear_site(lattice, props, site_new, site_one, prop_sum, active_indices);

            // new site can adsorb
            lattice->edges[site_new] = 'a';

        }

        lattice->edges[site_one] = 'd';

        return true;
    } // ADSORPTION 
    else if(reaction.type == Type::DESORPTION) {
        
        assert(lattice->edges.find(site_one) != lattice->edges.end());
        assert(lattice->sites[site_one].species == reaction.reactants[0]);
        assert(site_two == SITE_HOMOGENEOUS);

        lattice->sites[site_one].species = SPECIES_EMPTY;
        clear_site(lattice, props, site_one, std::optional<int> (), prop_sum, active_indices);

        lattice->edges[site_one] = 'a';
        
        if(is_add_sites) {
            // delete site
            lattice->delete_site(site_one);            
        }
  
        return true;

    } // DESORPTION
    else if((reaction.type == Type::HOMOGENEOUS_SOLID) || 
             (reaction.type == Type::REDUCTION) ||
             (reaction.type == Type::OXIDATION)){
        
        if(reaction.number_of_reactants == 1) {
            assert(lattice->sites[site_one].species == reaction.reactants[0]);
            assert(site_two == SITE_SELF_REACTION);
            lattice->sites[site_one].species = reaction.products[0];

            clear_site(lattice, props, site_one, std::optional<int> (), prop_sum, active_indices);

            if((lattice->edges.find(site_one) != lattice->edges.end()) && (lattice->edges[site_one] == 'd')) {
                
                clear_site_helper(props, site_one, SITE_HOMOGENEOUS, prop_sum, active_indices);
               
            }

        }
        else {
            // randomly assign products to sites 
            double r1 = sampler.generate();
            if(r1 <= 0.5) {
                lattice->sites[site_one].species = reaction.products[0];
                lattice->sites[site_two].species = reaction.products[1];
            }
            else {
                lattice->sites[site_one].species = reaction.products[1];
                lattice->sites[site_two].species = reaction.products[0];
            }

            if(lattice->edges.find(site_one) != lattice->edges.end()) {
                
                if((lattice->edges[site_one] == 'a') && (lattice->sites[site_one].species != SPECIES_EMPTY)) {
                    lattice->edges[site_one] = 'd';
                    clear_site_helper(props, site_one, SITE_HOMOGENEOUS, prop_sum, active_indices);
                }
                if(lattice->edges[site_one] == 'd') {
                    if(lattice->sites[site_one].species == SPECIES_EMPTY) {
                        lattice->edges[site_one] = 'a';

                        // remove site above from adsorb
                        if(is_add_sites) {
                            std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1};
                            int site_above = lattice->loc_map[key];

                            lattice->edges.erase(site_above);
                            clear_site_helper(props, site_above, SITE_HOMOGENEOUS, prop_sum, active_indices);
                        }
                    }
                    else {
                        clear_site_helper(props, site_one, SITE_HOMOGENEOUS, prop_sum, active_indices);
                    }
                }
               
            }

            if(lattice->edges.find(site_two) != lattice->edges.end()) {
                
                if((lattice->edges[site_two] == 'a') && (lattice->sites[site_two].species != SPECIES_EMPTY)) {
                    lattice->edges[site_two] = 'd';
                    clear_site_helper(props, site_two, SITE_HOMOGENEOUS, prop_sum, active_indices);
                }
                if(lattice->edges[site_two] == 'd') {
                    if(lattice->sites[site_two].species == SPECIES_EMPTY) {
                        lattice->edges[site_two] = 'a';

                        // remove site above from adsorb
                        if(is_add_sites) {
                            std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[site_two].i, lattice->sites[site_two].j, lattice->sites[site_two].k + 1};
                            int site_above = lattice->loc_map[key];

                            lattice->edges.erase(site_above);
                            clear_site_helper(props, site_above, SITE_HOMOGENEOUS, prop_sum, active_indices);
                        }
                    }
                    else {
                        clear_site_helper(props, site_two, SITE_HOMOGENEOUS, prop_sum, active_indices);
                    }
                }
               
            }


            clear_site(lattice, props, site_one, site_two, prop_sum, active_indices);
            clear_site(lattice, props, site_two, std::optional<int> (), prop_sum, active_indices);
            
        }

        return false;
    } // HOMOGENEOUS_SOLID OR REDUCTION OR OXIDATION 
    else if(reaction.type == Type::DIFFUSION) {

        int empty_site;
        int other_site;
        if(lattice->sites[site_one].species == SPECIES_EMPTY) {
            empty_site = site_one;
            other_site = site_two;
        } 
        else {
            assert(lattice->sites[site_two].species == SPECIES_EMPTY);
            empty_site = site_two;
            other_site = site_one;
        }
        lattice->sites[other_site].species = SPECIES_EMPTY;

        if(reaction.products[0] != SPECIES_EMPTY) {
            lattice->sites[empty_site].species = reaction.products[0];
        } 
        else {  
            lattice->sites[empty_site].species = reaction.products[1];
        }

        clear_site(lattice, props, site_one, site_two, prop_sum, active_indices);
        clear_site(lattice, props, site_two, std::optional<int> (), prop_sum, active_indices);

        if(lattice->edges.find(empty_site) != lattice->edges.end()) {
                
                if(lattice->edges[empty_site] == 'a'){
                   lattice->edges[empty_site] = 'd';
                    clear_site_helper(props, empty_site, SITE_HOMOGENEOUS, prop_sum, active_indices);
                }      
        }

        if(lattice->edges.find(other_site) != lattice->edges.end()) {
                
                if(lattice->edges[other_site] == 'd'){
                   lattice->edges[other_site] = 'a';

                    // remove site above from adsorb
                    if(is_add_sites) {
                        std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[other_site].i, lattice->sites[other_site].j, lattice->sites[other_site].k + 1};
                        int site_above = lattice->loc_map[key];

                        lattice->edges.erase(site_above);
                        clear_site_helper(props, site_above, SITE_HOMOGENEOUS, prop_sum, active_indices);
                   }
                   
                }      
        }       
        
        return false;
    } // DIFFUSION
    else {
        assert(false);
    } 
    
    
} // update_state() lattice

bool LatticeReactionNetwork::update_propensities(Lattice *lattice,
                        std::function<void(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props)> 
                        update_function, int next_reaction, int site_one, int site_two, 
                        std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props) {

    LatticeReaction reaction = reactions[next_reaction]; 

    if(reaction.type == Type::ADSORPTION) {

        // state already updated
        assert(lattice->sites[site_one].species == reaction.products[0]);
        assert(site_two == SITE_HOMOGENEOUS);

        relevant_react(lattice, update_function, site_one, std::optional<int> (), props);

        if(is_add_sites) {
            // also update new added site
            std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1};
            int site_new = lattice->loc_map[key];

            relevant_react(lattice, update_function, site_new, site_one, props);

        }
        
        return true;
    } // ADSORPTION 
    if(reaction.type == Type::DESORPTION) {
        
        // reaction already happended
        assert(lattice->sites[site_one].species == SPECIES_EMPTY);
        assert(site_two == SITE_HOMOGENEOUS);

        relevant_react(lattice, update_function, site_one, std::optional<int> (), props);

        if(is_add_sites) {
            // also update site above it that can no longer adsorb 
            std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1};
            int site_new = lattice->loc_map[key];

            relevant_react(lattice, update_function, site_new, site_one, props);

        }

        return true;
    } // DESORPTION
    else if((reaction.type == Type::HOMOGENEOUS_SOLID) || 
             (reaction.type == Type::REDUCTION) ||
             (reaction.type == Type::OXIDATION)) {
        
        if(reaction.number_of_reactants == 1) {
            assert(lattice->sites[site_one].species == reaction.products[0]);
            assert(site_two == SITE_SELF_REACTION);

            relevant_react(lattice, update_function, site_one, std::optional<int> (), props);
        }
        else {
    
            relevant_react(lattice, update_function, site_one, site_two, props);
            relevant_react(lattice, update_function, site_two, std::optional<int> (), props);
        }
        return false;
    } // HOMOGENEOUS_SOLID
    else if(reaction.type == Type::DIFFUSION) {

        relevant_react(lattice, update_function, site_one, site_two, props);
        relevant_react(lattice, update_function, site_two, std::optional<int> (), props);
        
        return false;
    } // DIFFUSION
    else {
        assert(false);
    }


} //update_propensities() lattice

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::clear_site(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int site, std::optional<int> ignore_neighbor, 
                        long double &prop_sum, int &active_indices) {

    
    assert(site != SITE_HOMOGENEOUS);
    // reset or initiate site combos
    for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++ ) {
        int neighbor = lattice->idneigh[site][neigh];

        if(!ignore_neighbor || (ignore_neighbor && neighbor != ignore_neighbor.value())) {
            clear_site_helper(props, site, neighbor, prop_sum, active_indices);
        }
    }

    // self reactions 
    clear_site_helper(props, site, SITE_SELF_REACTION, prop_sum, active_indices);
    

} // clear_site

/* ---------------------------------------------------------------------- */

// deal with active_indices
void LatticeReactionNetwork::clear_site_helper(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site_one, int site_two, long double &prop_sum,
                        int &active_indices) {

    std::string combo = make_string(site_one, site_two);
            
    // check if first time key has been added
    if(props.find(combo) == props.end()) {
        // key not found
        props[combo].reserve(reactions.size());
    }
    else {
        // already exists, clear vector to update
        std::vector<std::pair<double, int>> vec = props[combo];
        prop_sum -= sum_row(combo, props);

        if(prop_sum < 0) {
            assert(false);
        }

        active_indices -= props[combo].size();

        props[combo].clear();
    }

} // clear_site_helper

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::relevant_react(Lattice *lattice, std::function<void(LatticeUpdate lattice_update, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props)> update_function,
                                            int site, std::optional<int> ignore_neighbor, 
                                            std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props) {

    // all reactions related to central site 
    std::vector<int> potential_reactions = dependents[lattice->sites[site].species]; 

    // compute and add new propensities 
    for(size_t i = 0; i < potential_reactions.size(); i++ ) {

        unsigned long int reaction_id = potential_reactions[i];

        LatticeReaction reaction = reactions[reaction_id]; 

        // single reactant (must be of type lattice)
        if(reaction.number_of_reactants == 1 && reaction.phase_reactants[0] == Phase::LATTICE) {

            // if reaction produces a solid product make sure on edge 
            if(reaction.type == Type::DESORPTION) {
                //if(std::find(lattice->can_desorb.begin(), lattice->can_desorb.end(), site) != lattice->can_desorb.end()) {
                    if((lattice->edges.find(site) != lattice->edges.end()) &&
                    (lattice->edges[site] == 'd')) {
                    // in can_desorb vector
                    
                    double new_propensity = compute_propensity(1, 0, reaction_id, lattice);

                    update_function(LatticeUpdate {
                    .index = reaction_id,
                    .propensity = new_propensity,
                    .site_one = site,
                    .site_two = SITE_HOMOGENEOUS}, props); 
                }
            }
            else {
                // oxidation / reduction / homogeneous solid with one reactant 
                double new_propensity = compute_propensity(1, 0, reaction_id, lattice, site);

                update_function(LatticeUpdate {
                    .index = reaction_id,
                    .propensity = new_propensity,
                    .site_one = site,
                    .site_two = SITE_SELF_REACTION}, props); 
            }

        } // single reactant
        // two reactants (oxidation / reduction / homogeneous solid)
        else if(reaction.number_of_reactants == 2 && reaction.type != Type::ADSORPTION) {

            if(reaction.type == Type::OXIDATION && lattice->sites[site].species == 3) {
                assert(true);
            }
            
            int site_reactant_id; 
            int other_reactant_id; 
            
            if(lattice->sites[site].species == reaction.reactants[0]) {
                site_reactant_id = 0;
                other_reactant_id = 1;
            }
            else {
                site_reactant_id = 1;
                other_reactant_id = 0;
            }
            
            // check phase of site reactant 
            if(reaction.phase_reactants[site_reactant_id] == Phase::LATTICE) {

                // make sure neighbor is relevant 
                for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++) {
                    int neighbor = lattice->idneigh[site][neigh];
                    
                    if(!ignore_neighbor || (ignore_neighbor && neighbor != ignore_neighbor.value())) {
                        if(lattice->sites[neighbor].species == reaction.reactants[other_reactant_id]) {

                            double new_propensity = compute_propensity(1, 1, reaction_id, lattice);

                            update_function(LatticeUpdate {
                                            .index = reaction_id,
                                            .propensity = new_propensity,
                                            .site_one = site,
                                            .site_two = neighbor}, props); 
                        }
                    } // igore_neighbor

                } // for neigh
            }
            
        } // two reactants

    }

} // relevant_react

/* ---------------------------------------------------------------------- */

std::string LatticeReactionNetwork::make_string(int site_one, int site_two) {

    return (site_one > site_two) ? std::to_string(site_one) + "." + 
    std::to_string(site_two) : std::to_string(site_two) + "." + std::to_string(site_one);

} // make_string

std::string LatticeReactionNetwork::make_string(std::vector<int> vec) {
    std::sort(vec.begin(), vec.end());
    std::string hash = "";
    
    for(int i = 0; i < vec.size(); i++) {
        hash = hash + std::to_string(vec[i]);
    }
    
    return hash;

} // make_string

/* ---------------------------------------------------------------------- */

double LatticeReactionNetwork::sum_row(std::string hash, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props) {
    
    double sum = 0;    
    for(auto it = props[hash].begin(); it != props[hash].end(); it++) {
        sum += it->first;
    }
    return sum;
}

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::init_reaction_network(SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database, Lattice *lattice)

    {
    // collecting reaction network metadata
    SqlStatement<MetadataSql> metadata_statement (reaction_network_database);
    SqlReader<MetadataSql> metadata_reader (metadata_statement);

    std::optional<MetadataSql> maybe_metadata_row = metadata_reader.next();

    if (! maybe_metadata_row.has_value()) {
        std::cerr << time_stamp()
                  << "no metadata row\n";

        std::abort();
    }

    MetadataSql metadata_row = maybe_metadata_row.value();

    // setting reaction network factors
    SqlStatement<FactorsSql> factors_statement (initial_state_database);
    SqlReader<FactorsSql> factors_reader (factors_statement);

    // TODO: make sure this isn't nothing
    FactorsSql factors_row = factors_reader.next().value();
    factor_two = factors_row.factor_two;
    factor_duplicate = factors_row.factor_duplicate;

    // loading intial state
    initial_state.resize(metadata_row.number_of_species);

    SqlStatement<InitialStateSql> initial_state_statement (initial_state_database);
    SqlReader<InitialStateSql> initial_state_reader (initial_state_statement);

    int species_id;
    while(std::optional<InitialStateSql> maybe_initial_state_row =
          initial_state_reader.next()) {

        InitialStateSql initial_state_row = maybe_initial_state_row.value();
        species_id = initial_state_row.species_id;
        initial_state[species_id] = initial_state_row.count;
    }
    reactions.reserve(metadata_row.number_of_reactions);
    initial_propensities.resize(metadata_row.number_of_reactions);

    // loading reactions
    // vectors are default initialized to empty.
    // it is "cleaner" to resize the default vector than to
    // drop it and reinitialize a new vector.

    fill_reactions(reaction_network_database);
    assert(reactions.size() == metadata_row.number_of_reactions);

    // computing initial propensities
    for (unsigned long int i = 0; i < initial_propensities.size(); i++) {
        initial_propensities[i] = compute_propensity(initial_state, i, lattice);
    }

    std::cerr << time_stamp() << "computing dependency graph...\n";
    compute_dependents();
    std::cerr << time_stamp() << "finished computing dependency graph\n";

};

void LatticeReactionNetwork::fill_reactions(SqlConnection &reaction_network_database) {
    
    SqlStatement<LGMCReactionSql> reaction_statement (reaction_network_database);
    SqlReader<LGMCReactionSql> reaction_reader (reaction_statement);

    // reaction_id is lifted so we can do a sanity check after the
    // loop.  Make sure size of reactions vector, last reaction_id and
    // metadata number_of_reactions are all the same

    // read in Reactions
    while(std::optional<LGMCReactionSql> maybe_reaction_row = reaction_reader.next()) {

        LGMCReactionSql reaction_row = maybe_reaction_row.value();

        LatticeReaction reaction;

        reaction.number_of_reactants = static_cast<uint8_t> (reaction_row.number_of_reactants);
        reaction.number_of_products = static_cast<uint8_t> (reaction_row.number_of_products);
        reaction.reactants[0] = reaction_row.reactant_1;
        reaction.reactants[1] = reaction_row.reactant_2;
        reaction.products [0]= reaction_row.product_1;
        reaction.products[1] = reaction_row.product_2;
        reaction.rate = reaction_row.rate;
        reaction.dG = reaction_row.dG;
        reaction.prefactor = reaction_row.prefactor;
        reaction.reorganization_energy = reaction_row.reorganization_energy;
        reaction.electron_tunneling_coefficient = reaction_row.electron_tunneling_coefficient;
        reaction.charge_transfer_coefficient = reaction_row.charge_transfer_coefficient;

        if (strcmp(reinterpret_cast<const char *>(reaction_row.phase_reactant_1), "L") == 0) {
            reaction.phase_reactants[0] = Phase::LATTICE;
        } else if (strcmp(reinterpret_cast<const char *>(reaction_row.phase_reactant_1), "S") == 0) {
            reaction.phase_reactants[0] = Phase::SOLUTION;
        } else {
            reaction.phase_reactants[0] = Phase::NONE;
        }

        if (strcmp(reinterpret_cast<const char *>(reaction_row.phase_reactant_2), "L") == 0) {
            reaction.phase_reactants[1] = Phase::LATTICE;
        } else if (strcmp(reinterpret_cast<const char *>(reaction_row.phase_reactant_2), "S") == 0) {
            reaction.phase_reactants[1] = Phase::SOLUTION;
        } else {
            reaction.phase_reactants[1] = Phase::NONE;
        }

        if (strcmp(reinterpret_cast<const char *>(reaction_row.phase_product_1), "L") == 0) {
            reaction.phase_products[0] = Phase::LATTICE;
        } else if (strcmp(reinterpret_cast<const char *>(reaction_row.phase_product_1), "S") == 0) {
            reaction.phase_products[0] = Phase::SOLUTION;
        } else {
            reaction.phase_products[0] = Phase::NONE;
        }

        if (strcmp(reinterpret_cast<const char *>(reaction_row.phase_product_2), "L") == 0) {
            reaction.phase_products[1] = Phase::LATTICE;
        } else if (strcmp(reinterpret_cast<const char *>(reaction_row.phase_product_2), "S") == 0) {
            reaction.phase_products[1] = Phase::SOLUTION;
        } else {
            reaction.phase_products[1] = Phase::NONE;
        }
        
        if (strcmp(reinterpret_cast<const char *>(reaction_row.type), "O") == 0) {
            reaction.type = Type::OXIDATION;
        }
        else if (strcmp(reinterpret_cast<const char *>(reaction_row.type), "R") == 0) {
            reaction.type = Type::REDUCTION;
        }
        else if (strcmp(reinterpret_cast<const char *>(reaction_row.type), "F") == 0) {
            reaction.type = Type::DIFFUSION;
        }
        else if (strcmp(reinterpret_cast<const char *>(reaction_row.type), "L") == 0) {
            reaction.type =  Type::HOMOGENEOUS_SOLID;
        }
        else if (strcmp(reinterpret_cast<const char *>(reaction_row.type), "S") == 0) {
            reaction.type =  Type::HOMOGENEOUS_ELYTE;
        }
        else if (strcmp(reinterpret_cast<const char *>(reaction_row.type), "D") == 0) {
            reaction.type =  Type::DESORPTION;
        }
        else if (strcmp(reinterpret_cast<const char *>(reaction_row.type), "A") == 0) {
            reaction.type =  Type::ADSORPTION;
        }

        reactions.push_back(reaction);
    }
}



/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::compute_dependents() {
    // initializing dependency graph

    dependents.resize(initial_state.size());

    for ( unsigned int reaction_id = 0; reaction_id <  reactions.size(); reaction_id++ ) {
        LatticeReaction reaction = reactions[reaction_id]; 

        if(reaction.type != 'L') {
            for ( int i = 0; i < reaction.number_of_reactants; i++ ) {
                int reactant_id = reaction.reactants[i];
                dependents[reactant_id].push_back(reaction_id);
            }
        }

    }

};

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::update_state(std::vector<int> &state, int reaction_index) {
    
    LatticeReaction reaction = reactions[reaction_index];

    for (int m = 0;
         m < reaction.number_of_reactants;
         m++) {
        if(reaction.phase_reactants[m] == Phase::SOLUTION) {
            assert(state[reactions[reaction_index].reactants[m]] != 0);
            state[reactions[reaction_index].reactants[m]]--;
        }
        
    }

    for (int m = 0;
         m < reactions[reaction_index].number_of_products;
         m++) {
        if(reaction.phase_products[m] == Phase::SOLUTION) state[reactions[reaction_index].products[m]]++;
    }
};

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::update_propensities(std::function<void(Update update)> update_function,
                                                 std::vector<int> &state, int next_reaction, Lattice *lattice) {

    LatticeReaction reaction = reactions[next_reaction]; 

    std::vector<int> species_of_interest;
    species_of_interest.reserve(4);

    for ( int i = 0; i < reaction.number_of_reactants; i++ ) {
        // make sure reactant is in the solution
        if(reaction.phase_reactants[i] == Phase::SOLUTION) {
            int reactant_id = reaction.reactants[i];
            species_of_interest.push_back(reactant_id);
        }
        
    }

    for ( int j = 0; j < reaction.number_of_products; j++ ) {
        // make sure product is in the solution
        if(reaction.phase_products[j] == Phase::SOLUTION) {
            int product_id = reaction.products[j];
            species_of_interest.push_back(product_id);
        }
        
    }

    for ( int species_id : species_of_interest ) {
        for ( unsigned int reaction_index : dependents[species_id] ) {

            double new_propensity = compute_propensity(state, reaction_index, lattice);

            update_function(Update {
                    .index = reaction_index,
                    .propensity = new_propensity});
        }
    }

}

/* ---------------------------------------------------------------------- */

double LatticeReactionNetwork::compute_propensity(std::vector<int> &state, int reaction_index, Lattice *lattice) {
    LatticeReaction reaction = reactions[reaction_index];

    double p, k;

    for(int i = 0; i < reaction.number_of_reactants; i++) {
        if(reaction.phase_reactants[i] == Phase::LATTICE) {
            return 0;
        }
    }

    for(int i = 0; i < reaction.number_of_products; i++) {
        if(reaction.phase_products[i] == Phase::LATTICE) {
            return 0;
        }
    }

    if(reaction.type == Type::OXIDATION || reaction.type == Type::REDUCTION) {
         bool reduction_in = (reaction.type == Type::REDUCTION) ? true : false;

        if (charge_transfer_style == ChargeTransferStyle::MARCUS) {
            k = get_marcus_rate_coefficient(reaction.dG, reaction.prefactor, reaction.reorganization_energy,
                                            reaction.electron_tunneling_coefficient, g_e, lattice->get_maxz(), temperature,
                                            reduction_in);
        } else {
            k = get_butler_volmer_rate_coefficient(reaction.dG, reaction.prefactor, reaction.charge_transfer_coefficient,
                                                   reaction.electron_tunneling_coefficient, g_e, lattice->get_maxz(), temperature,
                                                   reduction_in);    
        }

    } else {
        k = reaction.rate;
    }
    
    // one reactant
    if (reaction.number_of_reactants == 1)
        p = state[reaction.reactants[0]] * k;


    // two reactants
    else {
        if (reaction.reactants[0] == reaction.reactants[1])
            p = factor_duplicate
                * factor_two
                * state[reaction.reactants[0]]
                * (state[reaction.reactants[0]] - 1)
                * k;

        else
            p = factor_two
                * state[reaction.reactants[0]]
                * state[reaction.reactants[1]]
                * k;
    }
    if(p < 0){
        assert(false);
    }
    if(p < 1e-6) {
        return 0;
    }
    else {
        return p;
    }

}

/* ---------------------------------------------------------------------- */

LatticeWriteTrajectoriesSql LatticeReactionNetwork::history_element_to_sql(
    int seed,
    LatticeTrajectoryHistoryElement history_element) {
    return LatticeWriteTrajectoriesSql {
        .seed = seed,
        .step = history_element.step,
        .time = history_element.time,
        .reaction_id = history_element.reaction_id,
        .site_1 = history_element.site_1,
        .site_2 = history_element.site_2
    };
}

/* ---------------------------------------------------------------------- */

LatticeWriteStateSql LatticeReactionNetwork::state_history_element_to_sql(
    int seed,
    LatticeStateHistoryElement state_history_element) {

    return LatticeWriteStateSql {
        .seed = seed,
        .species_id = state_history_element.species_id,
        .quantity = state_history_element.quantity,
        .i = state_history_element.i,
        .j = state_history_element.j,
        .k = state_history_element.k
    };
}

/* ---------------------------------------------------------------------- */

WriteCutoffSql LatticeReactionNetwork::cutoff_history_element_to_sql(
    int seed,
    CutoffHistoryElement cutoff_history_element) {
        return WriteCutoffSql {
            .seed = seed,
            .step = cutoff_history_element.step,
            .time = cutoff_history_element.time
        };
}

/* ---------------------------------------------------------------------- */

bool LatticeReactionNetwork::read_state(SqlReader<LatticeReadStateSql> state_reader, 
                    std::map<int, LatticeState> &temp_seed_state_map,
                LatticeReactionNetwork &lattice_reaction_network,
                SeedQueue &temp_seed_queue) {
    
    bool read_interrupt_states = false;
    Lattice *initial_lattice = lattice_reaction_network.initial_lattice;
    int initial_latconst = initial_lattice->latconst;

    // create a default lattice for each simulation
    while (std::optional<unsigned long int> maybe_seed =
        temp_seed_queue.get_seed()){
        unsigned long int seed = maybe_seed.value();

        // Each LatticeState must have its own lattice to point to
        Lattice *default_lattice = new Lattice(initial_latconst, 
        initial_lattice->xlo/initial_latconst, initial_lattice->xhi/initial_latconst, 
        initial_lattice->ylo/initial_latconst, initial_lattice->yhi/initial_latconst, 
        initial_lattice->zlo/initial_latconst, initial_lattice->zhi/initial_latconst,
        initial_lattice->is_xperiodic, initial_lattice->is_yperiodic, initial_lattice->is_zperiodic);

        LatticeState default_state = {lattice_reaction_network.initial_state, default_lattice};

        temp_seed_state_map.insert(std::make_pair(seed, default_state));

    }

    // TODO: add in for dynamic lattice
    while (std::optional<LatticeReadStateSql> maybe_state_row = state_reader.next()){
        read_interrupt_states = true;

        LatticeReadStateSql state_row = maybe_state_row.value();
        // determine if in lattice or homogeneous region
        if(state_row.i == SITE_HOMOGENEOUS) {
            // set the quantity of the homogeneous region vector associated with that seed
            temp_seed_state_map[state_row.seed].homogeneous[state_row.species_id] = state_row.quantity;
        }
        else {
            // update occupancy of that lattice site
            std::tuple<uint32_t, uint32_t, uint32_t> key = {state_row.i, state_row.j, state_row.k};
            int site_id = temp_seed_state_map[state_row.seed].lattice->loc_map[key];
            temp_seed_state_map[state_row.seed].lattice->sites[site_id].species = state_row.species_id;

            // check if can absorb, if so make sure edges now 'd'
            if(temp_seed_state_map[state_row.seed].lattice->sites[site_id].can_adsorb) {
                temp_seed_state_map[state_row.seed].lattice->edges[site_id] = 'd';
            }
        }
    }


    return read_interrupt_states;  

}

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::read_trajectories(SqlReader<LatticeReadTrajectoriesSql> trajectory_reader, 
                           std::map<int, LatticeState> &temp_seed_state_map, 
                           std::map<int, int> &temp_seed_step_map, 
                           std::map<int, double> &temp_seed_time_map,
                           LatticeReactionNetwork &lattice_reaction_network) {

                            std::cout << 'Inside read_trajectories' << std::endl;

}

void LatticeReactionNetwork::store_state_history(std::vector<LatticeStateHistoryElement> &state_packet,
    LatticeState &state, LatticeReactionNetwork &lattice_reaction_network, unsigned long int &seed) {

    // Lattice site
    for (auto site : state.lattice->sites) {
        
        state_packet.push_back(LatticeStateHistoryElement{
            .seed = seed,
            .species_id = site.second.species,
            .quantity = 1,
            .i = static_cast<int> (site.second.i),
            .j = static_cast<int> (site.second.j),
            .k = static_cast<int> (site.second.k)
        });
    }


    // Homogeneous region
    for (unsigned int i = 0; i < state.homogeneous.size(); i++) {
        state_packet.push_back(LatticeStateHistoryElement{
            .seed = seed,
            .species_id = static_cast<int>(i),
            .quantity = state.homogeneous[i],
            .i = SITE_HOMOGENEOUS,
            .j = SITE_HOMOGENEOUS,
            .k = SITE_HOMOGENEOUS
        });
    }


}

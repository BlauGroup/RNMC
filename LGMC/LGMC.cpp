#include "LGMC.h"
#include "memory.h"
#include "lattice.h"
#include "../core/dispatcher.h"
#include "LatSolver.h"

#include <string>
#include <fstream>
#include <iostream>
#include <getopt.h>


const int EMPTY_SITE = 0;
const int SELF_REACTION = -3;
const int GILLESPIE_SITE = -2;

const double TEMPERATURE = 298.15;     // In Kelvin
const double KB = 8.6173e-5;           // In eV/K
const double PLANCK = 4.1357e-15;      // In eV s

// NOT ACTUALLY CONSTANTS
// SHOULD PROBABLY BE COMPONENTS OF THE REACTION NETWORK PARAMETERS
const double TUNNEL_COEF = 1.2;        // Electron tunneling coefficient, in A^-1 
const double G_E = -2.15;              // Electron free energy, in eV


using namespace LGMC_NS;


double marcus_rate_coefficient(double base_dg, double reorganization_energy, double e_free, double distance, bool reduction) {

    double dg, dg_barrier, squared, kappa;

    if (reduction) {
        dg = base_dg - e_free;
    }
    else {
        dg = base_dg + e_free;
    }

    squared = 1 + dg / reorganization_energy;
    dg_barrier = reorganization_energy / 4 * squared * squared;
    kappa = std::exp(-1 * TUNNEL_COEF * distance);

    if (dg_barrier < 0) {
        return kappa * KB * TEMPERATURE / PLANCK;
    } else {
        return kappa * KB * TEMPERATURE / PLANCK * std::exp(-1 * dg_barrier / (KB * TEMPERATURE));
    }
}

LGMC::LGMC(SqlConnection &reaction_network_database, SqlConnection &initial_state_database, 
            LGMCParameters parameters) : sampler (Sampler(0)) {

    ReactionNetworkParameters react_parameters;
    LatticeReactionNetwork lattice_reaction_network = LatticeReactionNetwork(reaction_network_database, initial_state_database, react_parameters);
    react_net = &lattice_reaction_network;

    // create lattice
    lattice = new Lattice(parameters.latconst, parameters.boxxlo, parameters.boxxhi, parameters.boxylo,
                    parameters.boxyhi, parameters.boxzlo, parameters.boxzhi, parameters.xperiodic, parameters.yperiodic, parameters.zperiodic);

    initial_propensities = react_net->initial_propensities;
    initial_state = react_net->initial_state;
} // LGMC()

/* ---------------------------------------------------------------------- */

LGMC::~LGMC()
{
    // deal with pointers 

    delete lattice;
} // ~LGMC()


/* ---------------------------------------------------------------------- 
    Only calls this function if necessary reactants are on lattice sites
---------------------------------------------------------------------- */

double LGMC::compute_propensity(int site_one, int site_two, int num_one, int num_two, 
                                int react_id) {
    
    LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[react_id].get());

    double p;

    if(reaction->type == Type::OXIDATION || reaction->type == Type::REDUCTION) {
        double z;
        if (site_one < 0) {
            //TODO: FIX THIS. REALLY SHOULD BE THE LOWEST Z OF ANY SITE THAT CAN ADSORB
            z = lattice->zhi;
        }
        else {
            z = lattice->sites[site_one].z;
        }

        bool reduction = (reaction->type == Type::REDUCTION);

        p = num_one * marcus_rate_coefficient(reaction->dG, reaction->reorganization_energy, G_E, z, reduction);

    }

    // one reactant
    else if (reaction->number_of_reactants == 1)
        p = num_one * reaction->rate;

    // two reactants
    else {
        if (reaction->reactants[0] == reaction->reactants[1])
            p = react_net->factor_duplicate
                * react_net->factor_two
                * num_one
                * (num_one - 1)
                * reaction->rate;

        else
            p = react_net->factor_two
                * num_one
                * num_two
                * reaction->rate;
    }

    assert(p != 0);

    return p;
    
} // compute_propensity() 

/* ---------------------------------------------------------------------- */

bool LGMC::update_state(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int next_reaction, int site_one, int site_two, double &prop_sum) {

    LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[next_reaction].get());

    if(reaction->type == Type::ADSORPTION) {

        assert(lattice->sites[site_one].species == EMPTY_SITE);
        assert(lattice->sites[site_two].species == GILLESPIE_SITE);

        for(int i = 0; i < reaction->number_of_reactants; i++) {
            if(reaction->phase_reactants[i] == Phase::LATTICE) {

                // update site
                lattice->sites[site_one].species = reaction->products[0];
                clear_site(props, site_one, std::optional<int> (), prop_sum);
            }
        }

        return true;
    } // ADSORPTION 
    else if(reaction->type == Type::DESORPTION) {
        
        assert(lattice->sites[site_one].species == reaction->reactants[0]);
        assert(lattice->sites[site_two].species == SELF_REACTION);
        lattice->sites[site_one].species = EMPTY_SITE;

        clear_site(props, site_one, std::optional<int> (), prop_sum);

        return true;

    } // DESORPTION
    else if(reaction->type == Type::HOMOGENEOUS_SOLID) {
        
        if(reaction->number_of_reactants == 1) {
            assert(lattice->sites[site_one].species == reaction->reactants[0]);
            assert(lattice->sites[site_two].species == SELF_REACTION);
            lattice->sites[site_one].species = reaction->products[0];

            clear_site(props, site_one, std::optional<int> (), prop_sum);

        }
        else {
            // randomly assign products to sites 
            double r1 = sampler.generate();
            if(r1 <= 0.5) {
                lattice->sites[site_one].species = reaction->products[0];
                lattice->sites[site_two].species = reaction->products[1];
            }
            else {
                lattice->sites[site_one].species = reaction->products[1];
                lattice->sites[site_two].species = reaction->products[0];
            }

            clear_site(props, site_one, site_two, prop_sum);
            clear_site(props, site_two, std::optional<int> (), prop_sum);
            
        }

        return false;
    } // HOMOGENEOUS_SOLID
    else if(reaction->type == Type::DIFFUSION) {

        assert(lattice->sites[site_two].species == EMPTY_SITE);
        lattice->sites[site_one].species = EMPTY_SITE;

        if(reaction->products[0] != EMPTY_SITE) {
            lattice->sites[site_two].species = reaction->products[0];
        } 
        else {  
            lattice->sites[site_two].species = reaction->products[1];
        }

        clear_site(props, site_one, site_two, prop_sum);
        clear_site(props, site_two, std::optional<int> (), prop_sum);
        
        return false;
    } // DIFFUSION
    
    
} // update_state() lattice

void LGMC::update_propensities(std::vector<int> &state,
                        std::function<void(LatticeUpdate lattice_update)> 
                        update_function, int next_reaction, int site_one, int site_two) {

    LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[next_reaction].get());

    if(reaction->type == Type::ADSORPTION) {

        assert(lattice->sites[site_one].species == EMPTY_SITE);
        assert(lattice->sites[site_two].species == GILLESPIE_SITE);

        relevant_react(update_function, state, site_one, std::optional<int> ());

    } // ADSORPTION 
    if(reaction->type == Type::DESORPTION) {
        
        assert(lattice->sites[site_one].species == reaction->reactants[0]);
        assert(lattice->sites[site_two].species == SELF_REACTION);

        relevant_react(update_function, state, site_one, std::optional<int> ());


    } // DESORPTION
    else if(reaction->type == Type::HOMOGENEOUS_SOLID) {
        
        if(reaction->number_of_reactants == 1) {
            assert(lattice->sites[site_one].species == reaction->reactants[0]);
            assert(lattice->sites[site_two].species == SELF_REACTION);

            relevant_react(update_function, state, site_one, std::optional<int> ());
        }
        else {
    
            relevant_react(update_function, state, site_one, site_two);
            relevant_react(update_function, state, site_two, std::optional<int> ());
        }

    } // HOMOGENEOUS_SOLID
    else if(reaction->type == Type::DIFFUSION) {

        assert(lattice->sites[site_two].species == EMPTY_SITE);

        relevant_react(update_function, state, site_one, site_two);
        relevant_react(update_function, state, site_two, std::optional<int> ());
        
    } // DIFFUSION


} //update_propensities() lattice

/* ---------------------------------------------------------------------- */
/* ----------------*/
/* TODO: FOR EVAN */
/* ----------------*/
void LGMC::update_adsorption() {

    assert(false);

} // update_adsorption

/* ---------------------------------------------------------------------- */

void LGMC::clear_site(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int site, std::optional<int> ignore_neighbor, double &prop_sum) {

    // reset or initiate site combos
    for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++ ) {
        int neighbor = lattice->idneigh[site][neigh];

        if(!ignore_neighbor || (ignore_neighbor && neighbor != ignore_neighbor)) {
            clear_site_helper(props, site, neighbor, prop_sum);
        }
    }

    if(lattice->sites[site].can_adsorb) {
        // can interact with gillespie
        clear_site_helper(props, site, GILLESPIE_SITE, prop_sum);
    }

    // self reactions 
    if(site != EMPTY_SITE && site != GILLESPIE_SITE) {
        clear_site_helper(props, site, SELF_REACTION, prop_sum);
    }
    

} // clear_site

/* ---------------------------------------------------------------------- */

// deal with active_indicies
void LGMC::clear_site_helper(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site_one, int site_two, double &prop_sum) {

    std::string combo = make_string(site_one, site_two);
            
    // check if first time key has been added
    if(props.find(combo) == props.end()) {
        // key not found
        props[combo].reserve(react_net->dependents.size());
    }
    else {
        // already exists, clear vector to update
        prop_sum -= sum_row(combo, props);

        props[combo].clear();
    }

} // clear_site_helper

/* ---------------------------------------------------------------------- */

void LGMC::relevant_react(std::function<void(LatticeUpdate lattice_update)> update_function, 
                            std::vector<int> &state, int site, std::optional<int> ignore_neighbor) {

    // all reactions related to central site 
    std::vector<int> &potential_reactions = react_net->dependents[lattice->sites[site].species]; 

    // compute and add new propensities 
    for(unsigned long int reaction_id = 0; reaction_id < static_cast<int> (potential_reactions.size()); reaction_id++ ) {

        LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[reaction_id].get());

        // single reactant (must be of type lattice)
        if(reaction->number_of_reactants == 1 && reaction->phase_reactants[0] == Phase::LATTICE) {

            // if reaction produces a solid product make sure on edge 
            if(reaction->type == Type::DESORPTION) {
                if(lattice->sites[site].can_adsorb) {
                    
                    double new_propensity = compute_propensity(site, SELF_REACTION, 1, 0, reaction_id);

                    update_function(LatticeUpdate {
                    .index = reaction_id,
                    .propensity = new_propensity,
                    .site_one = site,
                    .site_two = SELF_REACTION}); 
                }
            }

        } // single reactant
        // two reactants 
        else if(reaction->number_of_reactants == 2) {
            
            int site_reactant_id; 
            int other_reactant_id; 
            
            if(lattice->sites[site].species == reaction->reactants[0]) {
                site_reactant_id = 0;
                other_reactant_id = 1;
            }
            else {
                site_reactant_id = 1;
                other_reactant_id = 0;
            }
            
            // check phase of site reactant 
            if(reaction->phase_reactants[site_reactant_id] == Phase::LATTICE) {
                
                if(reaction->type == Type::ADSORPTION) {
                    
                    double new_propensity = compute_propensity(site, GILLESPIE_SITE, 1, state[other_reactant_id], reaction_id);

                    update_function(LatticeUpdate {
                                    .index = reaction_id,
                                    .propensity = new_propensity,
                                    .site_one = site,
                                    .site_two = GILLESPIE_SITE});

                } // ADSORPTION 

                // make sure neighbor is relevant 
                for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++) {
                    int neighbor = lattice->idneigh[site][neigh];
                    
                    if(!ignore_neighbor || (ignore_neighbor && neighbor != ignore_neighbor)) {
                        if(lattice->sites[neighbor].species == reaction->reactants[other_reactant_id]) {

                            double new_propensity = compute_propensity(site, neighbor, 1, 1, reaction_id);

                            update_function(LatticeUpdate {
                                            .index = reaction_id,
                                            .propensity = new_propensity,
                                            .site_one = site,
                                            .site_two = neighbor}); 
                        }
                    } // igore_neighbor

                } // for neigh
            }
            
        } // two reactants

    }

} // relevant_react

/* ---------------------------------------------------------------------- */

std::string LGMC::make_string(int site_one, int site_two) {

    return (site_one < site_two) ? std::to_string(site_one) + "." + 
    std::to_string(site_two) : std::to_string(site_two) + "." + std::to_string(site_one);

} // make_string


/* -------------------------------- Gillespie Updates ----------------------------- */
void LGMC::update_state(std::vector<int> &state, int reaction_index) {
    react_net->update_state(state, reaction_index);
}

void LGMC::update_propensities(std::function<void(Update update)> update_function,
                                std::vector<int> &state, int next_reaction) {
    react_net->update_propensities(update_function, state, next_reaction);

}

/* ---------------------------------------------------------------------- */

double LGMC::sum_row(std::string hash, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props) {
    
    double sum = 0;    
    for(auto it = props[hash].begin(); it != props[hash].end(); it++) {
        sum += it->first;
    }
    return sum;
}

/* ---------------------------------------------------------------------- */

TrajectoriesSql LGMC::history_element_to_sql(
    int seed,
    HistoryElement history_element) {
    return TrajectoriesSql {
        .seed = seed,
        .step = history_element.step,
        .reaction_id = history_element.reaction_id,
        .time = history_element.time
    };
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

void print_usage() {
    
    std::cout << "Usage: specify the following in the input file\n"
              << "reaction_database\n"
              << "initial_state_database\n"
              << "number_of_simulations\n"
              << "base_seed\n"
              << "step_cutoff\n"
              << "lattice_parameters\n";
    
} // print_usage()

void print_usage_lattice() {

}

/* ---------------------------------------------------------------------- */

int main(int argc, char **argv) {
    if (argc != 8) {
        print_usage();
        exit(EXIT_FAILURE);
    }


    struct option long_options[] = {
        {"reaction_database", required_argument, NULL, 1},
        {"initial_state_database", required_argument, NULL, 2},
        {"number_of_simulations", required_argument, NULL, 3},
        {"base_seed", required_argument, NULL, 4},
        {"thread_count", required_argument, NULL, 5},
        {"step_cutoff", optional_argument, NULL, 6},
        {"time_cutoff", optional_argument, NULL, 7},
        {"lattice_parameters", required_argument, NULL, 8},
        {NULL, 0, NULL, 0}
        // last element of options array needs to be filled with zeros
    };

    int c;
    int option_index = 0;

    char *reaction_database = nullptr;
    char *initial_state_database = nullptr;
    int number_of_simulations = 0;
    int base_seed = 0;
    int thread_count = 0;
    std::string lattice_params_file; 

    Cutoff cutoff = {
        .bound =  { .step =  0 },
        .type_of_cutoff = step_termination
    };


    while ((c = getopt_long_only(
                argc, argv, "",
                long_options,
                &option_index)) != -1) {

        switch (c) {

        case 1:
            reaction_database = optarg;
            break;

        case 2:
            initial_state_database = optarg;
            break;

        case 3:
            number_of_simulations = atoi(optarg);
            break;

        case 4:
            base_seed = atoi(optarg);
            break;

        case 5:
            thread_count = atoi(optarg);
            break;

        case 6:
            cutoff.bound.step = atoi(optarg);
            cutoff.type_of_cutoff = step_termination;
            break;

        case 7:
            cutoff.bound.time = atof(optarg);
            cutoff.type_of_cutoff = time_termination;
            break;

        case 8:
            lattice_params_file = optarg;
            break;

        default:
            // if an unexpected argument is passed, exit
            print_usage();
            exit(EXIT_FAILURE);
            break;

        }

    }

    // read in LGMC parameters from file 
    std::ifstream fin;
    fin.open(lattice_params_file);

    if(!fin.is_open()) {
        std::cout << "Failed to open file: " << lattice_params_file << "\n";
        exit(EXIT_FAILURE);
    }

    float latconst;                               
    float boxxlo,boxxhi,boxylo,                   
          boxyhi,boxzlo,boxzhi;                       
    bool xperiodic, yperiodic, zperiodic;


    std::cin >> latconst >> boxxlo >> boxxhi >> boxylo >> boxyhi >> boxzlo >> boxzhi
    >> xperiodic >> yperiodic >> zperiodic;

    if(std::cin.fail()) {
        std::cout << "Incorrect file arguments.\n";
        exit(EXIT_FAILURE);
    }



    LGMCParameters parameters{.latconst = latconst, .boxxlo = boxxlo, .boxxhi = boxxhi, 
                              .boxylo = boxyhi, .boxzlo = boxzlo, .boxzhi = boxzhi, 
                              .xperiodic = xperiodic, .yperiodic = yperiodic, .zperiodic = zperiodic};                               

    Dispatcher<
        TreeSolver,
        LGMC,
        LGMCParameters,
        TrajectoriesSql
        >

        dispatcher (
        reaction_database,
        initial_state_database,
        number_of_simulations,
        base_seed,
        thread_count,
        cutoff,
        parameters, 
        true
        );

    dispatcher.run_dispatcher();
    exit(EXIT_SUCCESS);


}
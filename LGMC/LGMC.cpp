#include "LGMC.h"
#include "memory.h"
#include "lattice.h"
#include "../core/dispatcher.h"
#include "LatSolver.h"

#include <string>
#include <fstream>
#include <iostream>
#include <getopt.h>


int EMPTY_SITE = 0;
int SELF_REACTION = -3;
int GILLESPIE_SITE = -2;

using namespace LGMC_NS;

LGMC::LGMC(int argc, char **argv){

    float latconst;                               
    int boxxlo,boxxhi,boxylo,                   
    boxyhi,boxzlo,boxzhi;                       
    float xperiodic,yperiodic,zperiodic;
    

    std::string reaction_database;
    std::string initial_state_database;
    int number_of_simulations = 0;
    int base_seed = 0;
    int steps;
    
    sampler = Sampler(base_seed);
    Cutoff cutoff = {
        .bound =  { .step =  0 },
        .type_of_cutoff = step_termination
    };


    // Error if not only one file
    if(argc != 2) {
        std::cout << "Only input one file.\n";
        print_usage();
        exit(EXIT_FAILURE);
    }
    
    std::string file_in = argv[1];
    std::ifstream fin;
    fin.open(file_in);

    if(!fin.is_open()) {
        std::cout << "Failed to open file: " << file_in << "\n";
        exit(EXIT_FAILURE);
    }
    
    std::cin >> latconst >> boxxlo >> boxxhi >> boxylo >> boxyhi >> boxzlo >> boxzhi
    >> xperiodic >> yperiodic >> zperiodic >> reaction_database >> initial_state_database
    >> number_of_simulations >> base_seed >> steps;

    if(std::cin.fail()) {
        std::cout << "Incorrect file arguments.\n";
        exit(EXIT_FAILURE);
    }

    cutoff.bound.step = steps;

    // create lattice
    lattice = new Lattice(latconst, boxxlo, boxxhi, boxylo,
                          boxyhi, boxzlo, boxzhi, xperiodic, yperiodic, zperiodic);

    prop_sum = 0;

    ReactionNetworkParameters parameters;
   
} // LGMC()

LGMC::~LGMC()
{
    // deal with pointers 

    delete lattice;
} // ~LGMC()

void LGMC::print_usage() {
    
    std::cout << "Usage: specify the following in the input file\n"
              << "latconst, boxxlo, boxxhi, boxylo, boxyhi, boxzlo, boxzhi" 
              << "xperiodic (true/false), yperiodic (true/false), zperiodic (true/false)\n"
              << "reaction_database\n"
              << "initial_state_database\n"
              << "number_of_simulations\n"
              << "base_seed\n"
              << "step_cutoff\n";
    
} // print_usage()

/* ---------------------------------------------------------------------- */

void LGMC::run() {


} // run()

/* ---------------------------------------------------------------------- 
    Only calls this function if necessary reactants are on lattice sites
---------------------------------------------------------------------- */

void LGMC::compute_propensity(int site_one, int site_two, int num_one, int num_two, 
                              int react_id) {
    
    LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[next_reaction].get());

    double p;

    // TODO: make function for electrochemical reactions (for Evan)
    if(reaction->type == Type::ELECTROCHEMICAL) {
        assert(false);
    }

    // one reactant
    if (reaction->number_of_reactants == 1)
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

    // add exisiting propensity
    std::string site_combo = make_string(site_one, site_two);
    props[site_combo].push_back(std::make_pair(p, react_id));

    assert(p != 0);
    
} // compute_propensity() 

/* ---------------------------------------------------------------------- */

bool LGMC::update_state(int next_reaction, int site_one, int site_two) {

    LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[next_reaction].get());

    if(reaction->type == Type::ADSORPTION) {

        assert(sites[site_one].species == EMPTY_SITE);
        assert(sites[site_two].species == GILLESPIE_SITE);

        for(int i = 0; i < reaction->number_of_reactants; i++) {
            if(reaction->phase_reactants[i] == Phase::LATTICE) {

                // update site
                lattice->sites[empty_site].species = reaction->products[0];

            }
        }

        return true;
    } // ADSORPTION 
    else if(reaction->type == Type::DESORPTION) {
        
        assert(sites[site_one].species == reaction->reactants[0]);
        assert(sites[site_two].species == SELF_REACTION);
        lattice->sites[site_one].species = EMPTY_SITE;

        return true;

    } // DESORPTION
    else if(reaction->type == Type::HOMOGENEOUS_SOLID) {
        
        if(reaction->number_of_reactants == 1) {
            assert(sites[site_one].species == reaction->reactants[0]);
            assert(sites[site_two].species == SELF_REACTION);
            lattice->sites[site_one].species = reaction->products[0];

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
            
        }

        return false;
    } // HOMOGENEOUS_SOLID
    else if(reaction->type == Type::DIFFUSION) {

        assert(lattice->sites[site_two].species == EMPTY_SITE);
        lattice->sites[site_one] == EMPTY_SITE;

        if(reaction->products[0] != EMPTY_SITE) {
            lattice->sites[site_two] = reaction->products[0];
        } 
        else {  
            lattice->sites[site_two] = reaction->products[1];
        }
        
        return false;
    } // DIFFUSION
    
    
} // update_state() lattice

void LGMC::update_propensities(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        std::function<void(LatticeUpdate lattice_update)> update_function,
                        int next_reaction) {

    LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[next_reaction].get());

    if(reaction->type == Type::ADSORPTION) {

        assert(sites[site_one].species == EMPTY_SITE);
        assert(sites[site_two].species == GILLESPIE_SITE);

        clear_site(props, site_one, std::optional<int> ());
        relevant_react(update_function, site_one, std::optional<int> ());

    } // ADSORPTION 
    else if(reaction->type == Type::DESORPTION) {
        
        assert(sites[site_one].species == reaction->reactants[0]);
        assert(sites[site_two].species == SELF_REACTION);

        // clear, update site 
        clear_site(props, site_one, std::optional<int> ());
        relevant_react(update_function, site_one, std::optional<int> ());


    } // DESORPTION
    else if(reaction->type == Type::HOMOGENEOUS_SOLID) {
        
        if(reaction->number_of_reactants == 1) {
            assert(sites[site_one].species == reaction->reactants[0]);
            assert(sites[site_two].species == SELF_REACTION);

            // clear, update site 
            clear_site(props, site_one, std::optional<int> ());
            relevant_react(update_function, site_one, std::optional<int> ());
        }
        else {
            
            // clear, update site 
            clear_site(props, site_one, site_two);
            relevant_react(props, site_one, site_two);

            clear_site(props, site_two, std::optional<int> ());
            relevant_react(update_function, site_two, std::optional<int> ());
        }

    } // HOMOGENEOUS_SOLID
    else if(reaction->type == Type::DIFFUSION) {

        assert(lattice->sites[site_two].species == EMPTY_SITE);

        // clear, update site 
        clear_site(props, site_one, site_two);
        relevant_react(update_function, site_one, site_two);

        clear_site(props, site_two, std::optional<int> ());
        relevant_react(update_function, site_two, std::optional<int> ());
        
    } // DIFFUSION


} //update_propesnities() lattice

/* ---------------------------------------------------------------------- */

void LGMC::clear_site(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int site, std::optional<int> ignore_neighbor) {

    // reset or initiate site combos
    for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++ ) {
        int neighbor = lattice->idneigh[site][neigh];

        if(!ignore_neighbor || (ignore_neighbor && neighbor != ignore_neighbor)) {
            clear_site_helper(props, site, neighbor);
        }
    }

    if(lattice->is_on_edge(site)) {
        // can interact with gillespie
        clear_site_helper(props, site, GILLESPIE_SITE);
    }

    // self reactions 
    if(site != EMPTY_SITE && site != GILLESPIE_SITE) {
        clear_site_helper(props, site, SELF_REACTION);
    }
    

} // clear_site

/* ---------------------------------------------------------------------- */

// deal with active_indicies
void LGMC::clear_site_helper(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props
                        int site_one, int site_two) {

    std::string combo = make_string(site_one, site_two);
            
    // check if first time key has been added
    if(props.find(combo) == props.end()) {
        // key not found
        props[combo].reserve(lat_dependents.size());
    }
    else {
        // already exists, clear vector to update
        prop_sum -= sum_row(combo);

        // TODO: make sure capacity does not get changed
        props[combo].clear();
    }

} // clear_site_helper

/* ---------------------------------------------------------------------- */

void LGMC::relevant_react(std::function<void(LatticeUpdate lattice_update)> update_function,
                                             int site, std::optional<int> ignore_neighbor) {

    // all reactions related to central site 
    std::vector<int> &potential_reactions = react_net->dependents[lattice->sites[site].species]; 

    // compute and add new propensities 
    for(int reaction_id = 0; reaction_id < static_cast<int> (potential_reactions.size()); reaction_id++ ) {

        LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[reaction_id].get());

        // single reactant (must be of type lattice)
        if(reaction->number_of_reactants == 1 && reaction->phase_reactants[0] == Phase::LATTICE) {

            // if reaction produces a solid product make sure on edge 
            if(reaction->type == Type::DESORPTION) {
                if(lattice->is_on_edge(site)) {
                    
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
                // handling asorption in another function, other reactant must be on lattice
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

int main(int argc, char **argv) { 

    Dispatcher<
        LatSolver,
        LGMC,
        ReactionNetworkParameters,
        TrajectoriesSql
        >

        dispatcher (
        reaction_database,
        initial_state_database,
        number_of_simulations,
        base_seed,
        0,
        cutoff,
        parameters, 
        false
        );

    dispatcher.run_dispatcher();


    LGMC *battery = new LGMC(argc, argv);
    battery->run();

    exit(EXIT_SUCCESS);

}
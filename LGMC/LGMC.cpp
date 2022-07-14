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

    Dispatcher<
        LatSolver,
        LatticeReactionNetwork,
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

void LGMC::update_propensity(int site_one, int site_two, LatticeReaction *reaction, int react_id, int num_one, int num_two) {
    
    // TODO: make function for electrochemical reactions (for Evan)
    if(reaction->type == Type::ELECTROCHEMICAL) {
        update_electrochemical();
    }

    double p;

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

    // add or change existing propensity 
    std::string site_combo = make_string(site_one, site_two);
    props[site_combo].push_back(std::make_pair(p, react_id));

    // update running sum 
    prop_sum += p;


    
} // update_propensity()

/* ---------------------------------------------------------------------- */

void LGMC::update(std::optional<int> site_one, std::optional<int> site_two, int reaction_id) {

    LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[reaction_id].get());

    if(site_one) {
        // site_one specified must be a lattice reaction 
        if(reaction->number_of_reactants == 1) {
            update_one_lattice_site(*site_one, reaction);

        } // one reactant
        else if(reaction->number_of_reactants == 2) {
            int sites_to_update = 0; 
            for(int i = 0; i < reaction->number_of_reactants; i++) {
                if(reaction->phase_reactants[i] == Phase::LATTICE) {
                    sites_to_update++;
                }
            }

            if(sites_to_update == 1) {
                update_one_lattice_site(*site_one, reaction);
            }
            else {
                update_two_lattice_sites(*site_one, *site_two, reaction);
            }

        } // two reactants 
    }
    else {/*
        // Gillespie reaction 
        for(int i = 0; i < reaction->number_of_reactants; i++) {
            if(reaction->phase_reactants[i] == Phase::SOLUTION) {
               react_net->decrease_species(reaction->reactants[i]);
            }
            else {
                //find site to become empty 
                int site_id = *react_net->lattice_sites[reaction->reactants[i]].end(); 
                react_net->lattice_sites[reaction->reactants[i]].pop_back();

                lattice->sites[site_id].species = EMPTY_SITE;
                update_site(site_id, std::optional<int> ());

            }
        } // reactants 

        for(int i = 0; i < reaction->num_of_reactants; i++) {
            if(reaction->phase_reactants[i] == Phase::SOLUTION) {
               react_new->increase_species(reaction->reactants[i]);
            }
            else {
                int site_id = *react_net->lattice_sites[0].end(); 
                react_net->lattice_sites[0].pop_back();
                lattice->sites[site_id].species = reaction->products[i];
                update_site(site_id, std::optional<int> ());

            }
            // update propensities 
        react_net->update_propensities(std::function<void(Update update)> update_function,
                             int next_reaction);
                             */
        } // products
        

    }

    
} // update()
/* ---------------------------------------------------------------------- */

void LGMC::update_one_lattice_site(int site, LatticeReaction *reaction) {
    bool update_reaction = false;
    // if involves solution reactant update that 
    if(reaction->number_of_reactants == 2) {
        for(int i = 0; i < 2; i++) {
            if(reaction->phase_reactants[i] == Phase::SOLUTION) {
                react_net->decrease_species(reaction->reactants[i]);
                update_reaction = true;
            }
        }
    }
    // update lattice site species depending on reaction products
    // TODO: improve double setting of lattice to empty if both solution products
    for(int i = 0; i < reaction->number_of_products; i++) {
        if(reaction->phase_products[i] == Phase::LATTICE) {
            lattice->sites[site].species = reaction->products[i];
        }
        else {
            // solution product 
            lattice->sites[site].species = EMPTY_SITE;
            react_net->increase_species(reaction->products[i]);
            update_reaction = true;
        }
    }
    if(update_reaction) {
        react_net->update_propensities(std::function<void(Update update)> update_function,
                        int next_reaction);
    }
    
    clear_site(site, std::optional<int> ());
    update_site(site, std::optional<int> ());

} // update_one_lattice_site()

/* ---------------------------------------------------------------------- */

void LGMC::update_two_lattice_sites(int site_one, int site_two, LatticeReaction *reaction) {

    // figure out how many solution products
    int solid_products = 0;
    for(int i = 0; i < reaction->number_of_products; i++) {
        if(reaction->phase_products[i] == Phase::SOLUTION) {
            solid_products++;
        }
    }

    if(solid_products == 2) {

        // set both sites to empty
        lattice->sites[site_one].species = EMPTY_SITE;
        lattice->sites[site_two].species = EMPTY_SITE;

        //update reaction network 
        for(int i = 0; i < reaction->number_of_products; i++) { 
            react_net->increase_species(reaction->products[i]);
        }

        react_net->update_propensities(std::function<void(Update update)> update_function,
                        int next_reaction);
    } // two solution products
    else if(solid_products == 1) {

        for(int i = 0; i < reaction->number_of_products; i++) {
            if(reaction->phase_products[i] == Phase::SOLUTION) {
                react_net->increase_species(reaction->products[i]);
            }
            else {
                // use random number to determine which site gets the product, other becomes empty
                double r1 = sampler.generate();
                if(r1 <= 0.5) {
                    lattice->sites[site_one].species = reaction->products[i];
                    lattice->sites[site_two].species = EMPTY_SITE;
                }
                else {
                    lattice->sites[site_two].species = reaction->products[i];
                    lattice->sites[site_one].species = EMPTY_SITE;
                }
            }
        }

        react_net->update_propensities(std::function<void(Update update)> update_function,
                        int next_reaction);

    } // one solution product
    else {
       // if diffusion, make sure empty site gets product 
       if(reaction->type == Type::DIFFUSION) {
            int empty_site = 0;
            int non_empty = 0;
            if(lattice->sites[site_one].species == EMPTY_SITE) {
                empty_site = site_one;
                non_empty = site_two;
            }
            else {
                empty_site = site_two;
                non_empty = site_two;
            }
            if(reaction->products[0] != EMPTY_SITE) {
                lattice->sites[empty_site].species = reaction->products[0];
            }
            else {
                lattice->sites[empty_site].species = reaction->products[1];
            }
            lattice->sites[non_empty].species = EMPTY_SITE;
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
    }  // no solution product 

    // clear both sites, have first site ignore its neighbor
    clear_site(site_one, site_two);
    clear_site(site_two, std::optional<int> ());

    update_site(site_one, site_two);
    update_site(site_two, std::optional<int> ());

} // update_two_lattice_sites()

/* ---------------------------------------------------------------------- */

void LGMC::clear_site(int site, std::optional<int> ignore_neighbor) {

    // reset or initiate site combos
    for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++ ) {
        int neighbor = lattice->idneigh[site][neigh];

        if(!ignore_neighbor || (ignore_neighbor && neighbor != ignore_neighbor)) {
            clear_site_helper(site, neighbor);
        }
    }

    if(lattice->is_on_edge(site)) {
        // can interact with gillespie
        clear_site_helper(site, GILLESPIE_SITE);
    }

    // self reactions 
    if(site != EMPTY_SITE && site != GILLESPIE_SITE) {
        clear_site_helper(site, SELF_REACTION);
    }
    

} // clear_site

/* ---------------------------------------------------------------------- */

void LGMC::clear_site_helper(int site_one, int site_two) {

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

void LGMC::relevant_react(int site) {

    // all reactions related to central site 
    std::vector<int> &potential_reactions = lat_dependents[lattice->sites[site].species]; 

    // compute and add new propensities 
    for(int reaction_id = 0; reaction_id < static_cast<int> (potential_reactions.size()); reaction_id++ ) {

        LatticeReaction *reaction = static_cast<LatticeReaction *> (react_net->reactions[reaction_id].get());


        // single reactant (must be of type lattice)
        if(reaction->number_of_reactants == 1 && reaction->phase_reactants[0] == Phase::LATTICE) {

            // if reaction produces a solid product make sure on edge 
            if(reaction->type == Type::ADSORPTION) {
                if(lattice->is_on_edge(site)) {
                    update_propensity(site, SELF_REACTION, reaction, reaction_id, 1, 0);
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
            if(reaction->phase_reactants[site_reactant_id] == Phase::SOLUTION) {
                break;
            }
            if(reaction->phase_reactants[other_reactant_id] == Phase::SOLUTION) {
                // make sure site on edge then you can proceed 
                if(lattice->is_on_edge(site)) update_propensity(site, GILLESPIE_SITE, reaction, reaction_id, 1, react_net->state[reaction->reactants[other_reactant_id]]);
            }
            else {
                
                // make sure neighbor is relevant 
                for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++) {
                    int neighbor = lattice->idneigh[site][neigh];
                        
                    if(lattice->sites[neighbor].species == reaction->reactants[other_reactant_id]) {
                        // if adsoprtion make sure on edge 
                        if(reaction->type == Type::ADSORPTION && !(lattice->is_on_edge(site) || lattice->is_on_edge(neighbor))) {
                            break;
                        }
                        update_propensity(site, neighbor, reaction, reaction_id, 1, 1);
                    }

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

/* ---------------------------------------------------------------------- */

int LGMC::sum_row(std::string hash) {
    
    int sum = 0;    
    for(auto it = props[hash].begin(); it != props[hash].end(); it++) {
        sum += it->first;
    }
    return sum;
}

/* ---------------------------------------------------------------------- */

int main(int argc, char **argv) { 

    LGMC *battery = new LGMC(argc, argv);
    battery->run();

    exit(EXIT_SUCCESS);

}
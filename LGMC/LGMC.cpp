#include "spktyp.h"
#include "LGMC.h"
#include "memory.h"
#include "lattice.h"
#include <string>
#include <fstream>

//#include "../GMC/GMC.cpp"
#include <iostream>
#include <getopt.h>


// types of possible reactions
enum {ELECTROCHEMICAL, DIFFUSION, CHEMICAL, DESORPTION};

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

    // TODO: fix up dispatcher 
    Dispatcher< TreeSolver, ReactionNetwork,
        ReactionNetworkParameters,
        TrajectoriesSql >

    dispatcher (reaction_database, initial_state_database,
        number_of_simulations, base_seed, 0,
        cutoff, parameters);
   
} // LGMC()

LGMC::~LGMC()
{
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

void LGMC::update_propensity(int site_one, int site_two, Reaction &reaction, int react_id) {
    
    // TODO: for interactions with electrolyte make sure possible (take into account distance??)
    // TODO: with reactions with electrolyte take in number 
    // TODO: make function for electrochemical reactions 

    double p;
    // one reactant
    if (reaction.number_of_reactants == 1) {
        p = reaction.rate;
    }

    // two reactants
    else {
       if (reaction.reactants[0] == reaction.reactants[1]) {
            p = react_net->factor_duplicate
                * react_net->factor_two
                * 2 * reaction.rate;
       }
            

        else {
            p = react_net->factor_two * reaction.rate;
        }
    }

    // add or change existing propensity 
    std::string site_combo = make_string(site_one, site_two);
    props[site_combo].push_back(std::make_pair(p, react_id));

    // update running sum 
    prop_sum += p;


    
} // update_propensity()

/* ---------------------------------------------------------------------- */

void LGMC::update(int site_one, int site_two) {

    // TODO: overcounting 

    if(site_one > -1) {
        relevant_react(site_one);       
    }
    if(site_two > -1) {
        relevant_react(site_two);
    }
    
}

/* ---------------------------------------------------------------------- */

void LGMC::clear_site(int site) {

    // reset or initiate site combos
    for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++ ) {
        int neighbor = lattice->idneigh[site][neigh];

        std::string combo = make_string(site, neighbor);
        
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
    }

    // reset for gillespie (-2) and empty site (-1)
    for(int i = -2; i < 0; i++) {
        std::string combo = std::to_string(site) + "." + std::to_string(i);

        if(props.find(combo) == props.end()) {
         props[combo].resize(lat_dependents.size());
        }
        else {
            prop_sum -= sum_row(combo);
            props[combo].clear();
        }

    }

} // clear_site

/* ---------------------------------------------------------------------- */

void LGMC::relevant_react(int site) {

    // TODO: add in 

    clear_site(site);

    // all reactions related to central site 
    std::vector<int> &potential_reactions = lat_dependents[lattice->sites[site].species]; 

    // compute and add new propensities 
    for(int reaction_id = 0; reaction_id < static_cast<int> (potential_reactions.size()); reaction_id++ ) {

        Reaction &reaction = react_net->reactions[reaction_id];

        // reaction with gillepsie 
        if(reaction.type == DESORPTION) {
            update_propensity(site, -2, reaction, reaction_id);
        }

        // single reaction 
        if(reaction.number_of_reactants == 1) {
            
            if(reaction.number_of_products == 1) {
                update_propensity(site, -1, reaction, reaction_id);
            }
            else {
                // two products make sure site is empty
                for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++) {
                    
                    int neighbor = lattice->idneigh[site][neigh];
                    
                    if(lattice->sites[neighbor].species == -1) {
                        // empty 
                        update_propensity(site, neighbor, reaction, reaction_id);
                    }


                } // for neigh
            }

        } // single reactant

        if(reaction.number_of_reactants == 2) {
    
            int other_reactant = (lattice->sites[site].species == reaction.reactants[0]) ? reaction.reactants[1] : reaction.reactants[0];
            // make sure neighbor is relevant 
            for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++) {
                int neighbor = lattice->idneigh[site][neigh];
                    
                if(lattice->sites[neighbor].species == other_reactant) {
                    update_propensity(site, neighbor, reaction, reaction_id);
                }

            } // for neigh
            
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
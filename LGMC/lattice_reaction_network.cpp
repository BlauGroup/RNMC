#include "lattice_reaction_network.h"
#include "memory.h"
#include "lattice.h"

#include <string>
#include <assert.h>


const int SPECIES_EMPTY = 0;
const int SITE_SELF_REACTION = -3;
const int SITE_GILLESPIE = -2;

const double TEMPERATURE = 298.15;     // In Kelvin
const double KB = 8.6173e-5;           // In eV/K
const double PLANCK = 4.1357e-15;      // In eV s

// NOT ACTUALLY CONSTANTS
// SHOULD PROBABLY BE COMPONENTS OF THE REACTION NETWORK PARAMETERS
const double TUNNEL_COEF = 1.2;        // Electron tunneling coefficient, in A^-1 
const double G_E = -2.15;              // Electron free energy, in eV


using namespace LGMC_NS;


double get_marcus_rate_coefficient(double base_dg, double reorganization_energy, double e_free, double distance, bool reduction) {

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

LatticeReactionNetwork::LatticeReactionNetwork(SqlConnection &reaction_network_database, SqlConnection &initial_state_database, 
            LGMCParameters parameters) : sampler (Sampler(0)) {

    // create lattice
    initial_lattice = new Lattice(parameters.latconst, parameters.boxxlo, parameters.boxxhi, parameters.boxylo,
                    parameters.boxyhi, parameters.boxzlo, parameters.boxzhi, parameters.xperiodic, parameters.yperiodic, parameters.zperiodic);

    init_reaction_network(reaction_network_database, initial_state_database, initial_lattice);

    is_add_sites = parameters.is_add_sites;
    temperature = parameters.temperature;
    g_e = parameters.g_e;

} // LatticeReactionNetwork()

/* ---------------------------------------------------------------------- */

LatticeReactionNetwork::~LatticeReactionNetwork()
{
    // deal with pointers 

    delete initial_lattice;
} // ~LatticeReactionNetwork()

/* ---------------------------------------------------------------------- 
    Global Updates
---------------------------------------------------------------------- */

void LatticeReactionNetwork::update_state(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        std::vector<int> &state, int next_reaction, 
                        std::optional<int> site_one, std::optional<int> site_two, 
                        double &prop_sum, int &active_indices) {
    
    if(site_one) {
        // update lattice state
        bool update_gillepsie = update_state(lattice, props, next_reaction, site_one.value(), 
                                            site_two.value(), prop_sum, active_indices);
        
        if(update_gillepsie) {
            update_state(state, next_reaction);
        }
    }
    
    else {
        // gillespie event happens
        update_state(state, next_reaction);

    }

    update_adsorp_state(lattice, props, prop_sum, active_indices);

} // update_state

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::update_propensities(Lattice *lattice, std::vector<int> &state, std::function<void(Update update)> update_function, 
                                        std::function<void(LatticeUpdate lattice_update)> lattice_update_function, 
                                        int next_reaction, std::optional<int> site_one, std::optional<int> site_two) {

    if(site_one) {
        // update lattice state
        bool update_gillepsie = update_propensities(lattice, lattice_update_function, next_reaction, site_one.value(), site_two.value());
        
        if(update_gillepsie) {
            update_propensities(update_function, state, next_reaction, lattice);
        }
    }
    
    else {
        // gillespie event happens
        update_propensities(update_function, state, next_reaction, lattice);

    }

    update_adsorp_props(lattice, lattice_update_function, state);

}

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::update_adsorp_state(Lattice *lattice, std::unordered_map<std::string, std::vector< std::pair<double, int> > > &props,
                                double &prop_sum, int &active_indicies) {
    // update only sites on the edge 
    for(size_t i = 0; i < lattice->can_adsorb.size(); i++) {
        int site = lattice->can_adsorb[i];
        clear_site_helper(props, site, SITE_GILLESPIE, prop_sum, active_indicies);
    }

} // update_adsorp_state

void LatticeReactionNetwork::update_adsorp_props(Lattice *lattice, std::function<void(LatticeUpdate lattice_update)> lattice_update_function, std::vector<int> &state) {

        // update only sites on the edge 
    for(size_t i = 0; i < lattice->can_adsorb.size(); i++) {
        int site = lattice->can_adsorb[i];

        // find relevant adsorption reactions
        std::vector<int> &potential_reactions = dependents[lattice->sites[site].species]; 

        for(unsigned long int reaction_id = 0; reaction_id < static_cast<int> (potential_reactions.size()); reaction_id++ ) {

            LatticeReaction reaction = reactions[reaction_id]; 

            if(reaction.type == Type::ADSORPTION) {

                int other_reactant_id = (reaction.reactants[0] == lattice->sites[site].species) ? 1 : 0;
                
                double new_propensity = compute_propensity(1, state[other_reactant_id], reaction_id, lattice);

                lattice_update_function(LatticeUpdate {
                    .index = reaction_id,
                    .propensity = new_propensity,
                    .site_one = site,
                    .site_two = SITE_GILLESPIE}); 
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

    double p;

    if(reaction.type == Type::REDUCTION || reaction.type == Type::OXIDATION) {
        bool reduction_in = (reaction.type == Type::REDUCTION) ? true : false;

        // TODO: make general for all types of periodicity (distance = .z)
        get_marcus_rate_coefficient(reaction.dG, reaction.reorganization, g_e, lattice->sites[site_id].z, reduction_in);
    }

    // one reactant
    if (reaction.number_of_reactants == 1)
        p = num_one * reaction.rate;


    // two reactants
    else {
        if (reaction.reactants[0] == reaction.reactants[1])
            p = factor_duplicate
                * factor_two
                * num_one
                * (num_one - 1)
                * reaction.rate;

        else
            p = factor_two
                * num_one
                * num_two
                * reaction.rate;
    }

    assert(p != 0);

    return p;
    
} // compute_propensity() 

/* ---------------------------------------------------------------------- */

bool LatticeReactionNetwork::update_state(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int next_reaction, int site_one, int site_two, 
                        double &prop_sum, int &active_indices) {

    LatticeReaction reaction = reactions[next_reaction]; 

    if(reaction.type == Type::ADSORPTION) {

        assert(lattice->sites[site_one].species == SPECIES_EMPTY);
        assert(site_two == SITE_GILLESPIE);
        assert(std::find(lattice->can_adsorb.begin(), lattice->can_adsorb.end(), site_one) != lattice->can_adsorb.end());

        // update site
        lattice->sites[site_one].species = reaction.products[0];
        clear_site(lattice, props, site_one, std::optional<int> (), prop_sum, active_indices);

        if(is_add_sites) {
            // TODO: make general for all types of periodicity (add_site coordinates)
            lattice->add_site(lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1,
            lattice->sites[site_one].x, lattice->sites[site_one].y, lattice->sites[site_one].z + lattice->get_latconst(), 
            true, true, true);
            
            std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1};
            int site_new = lattice->loc_map[key];

            clear_site(lattice, props, site_new, site_one, prop_sum, active_indices);

            // remove initial site from can_adsorb 
            std::erase(lattice->can_adsorb, site_one);

            // add initial_site to can_desorb
            lattice->can_desorb.push_back(site_one);

        }

        return true;
    } // ADSORPTION 
    else if(reaction.type == Type::DESORPTION) {
        
        assert(lattice->sites[site_one].species == reaction.reactants[0]);
        assert(site_two == SITE_SELF_REACTION);

        lattice->sites[site_one].species = SPECIES_EMPTY;
        clear_site(lattice, props, site_one, std::optional<int> (), prop_sum, active_indices);

        if(is_add_sites) {
            assert(std::find(lattice->can_desorb.begin(), lattice->can_desorb.end(), site_one) != lattice->can_desorb.end());

            // remove site from desorb vector 
            auto it = std::find(lattice->can_desorb.begin(), lattice->can_desorb.end(), site_one);
            std::erase(lattice->can_desorb, site_one);

            // add to adsorb vector 
            lattice->can_adsorb.push_back(site_one);

            // remove site above from adsorb vector
            std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1};
            int site_above = lattice->loc_map[key];

            // remove site above from adsorb vector
            std::erase(lattice->can_adsorb, site_above);
        }
  
        return true;

    } // DESORPTION
    else if(reaction.type == Type::HOMOGENEOUS_SOLID) {
        
        if(reaction.number_of_reactants == 1) {
            assert(lattice->sites[site_one].species == reaction.reactants[0]);
            assert(lattice->sites[site_two].species == SITE_SELF_REACTION);
            lattice->sites[site_one].species = reaction.products[0];

            clear_site(lattice, props, site_one, std::optional<int> (), prop_sum, active_indices);

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

            clear_site(lattice, props, site_one, site_two, prop_sum, active_indices);
            clear_site(lattice, props, site_two, std::optional<int> (), prop_sum, active_indices);
            
        }

        return false;
    } // HOMOGENEOUS_SOLID
    else if(reaction.type == Type::DIFFUSION) {

        assert(lattice->sites[site_two].species == SPECIES_EMPTY);
        lattice->sites[site_one].species = SPECIES_EMPTY;

        if(reaction.products[0] != SPECIES_EMPTY) {
            lattice->sites[site_two].species = reaction.products[0];
        } 
        else {  
            lattice->sites[site_two].species = reaction.products[1];
        }

        clear_site(lattice, props, site_one, site_two, prop_sum, active_indices);
        clear_site(lattice, props, site_two, std::optional<int> (), prop_sum, active_indices);
        
        return false;
    } // DIFFUSION
    else if(reaction.type == Type::REDUCTION || reaction.type == Type::OXIDATION) {
        assert(lattice->sites[site_one].species == reaction.reactants[0]);
        assert(lattice->sites[site_two].species == SPECIES_EMPTY);

        lattice->sites[site_one].species = reaction.products[0];
        clear_site(lattice, props, site_one, std::optional<int> (), prop_sum, active_indices);
    } // reduction or oxidation 
    
    
} // update_state() lattice

bool LatticeReactionNetwork::update_propensities(Lattice *lattice,
                        std::function<void(LatticeUpdate lattice_update)> 
                        update_function, int next_reaction, int site_one, int site_two) {

    LatticeReaction reaction = reactions[next_reaction]; 

    if(reaction.type == Type::ADSORPTION) {

        assert(lattice->sites[site_one].species == SPECIES_EMPTY);
        assert(site_two == SITE_GILLESPIE);

        relevant_react(lattice, update_function, site_one, std::optional<int> ());

        if(is_add_sites) {
            // also update new added site
            std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1};
            int site_new = lattice->loc_map[key];

            relevant_react(lattice, update_function, site_new, site_one);

        }
        
        return true;
    } // ADSORPTION 
    if(reaction.type == Type::DESORPTION) {
        
        assert(lattice->sites[site_one].species == reaction.reactants[0]);
        assert(site_two == SITE_SELF_REACTION);

        relevant_react(lattice, update_function, site_one, std::optional<int> ());

        if(is_add_sites) {
            // also update site above it that can no longer adsorb 
            std::tuple<uint32_t, uint32_t, uint32_t> key = {lattice->sites[site_one].i, lattice->sites[site_one].j, lattice->sites[site_one].k + 1};
            int site_new = lattice->loc_map[key];

            relevant_react(lattice, update_function, site_new, site_one);

        }

        return true;
    } // DESORPTION
    else if(reaction.type == Type::HOMOGENEOUS_SOLID) {
        
        if(reaction.number_of_reactants == 1) {
            assert(lattice->sites[site_one].species == reaction.reactants[0]);
            assert(site_two == SITE_SELF_REACTION);

            relevant_react(lattice, update_function, site_one, std::optional<int> ());
        }
        else {
    
            relevant_react(lattice, update_function, site_one, site_two);
            relevant_react(lattice, update_function, site_two, std::optional<int> ());
        }
        return false;
    } // HOMOGENEOUS_SOLID
    else if(reaction.type == Type::DIFFUSION) {

        assert(lattice->sites[site_two].species == SPECIES_EMPTY);

        relevant_react(lattice, update_function, site_one, site_two);
        relevant_react(lattice, update_function, site_two, std::optional<int> ());
        
        return false;
    } // DIFFUSION
    else if(reaction.type == Type::REDUCTION || reaction.type == Type::OXIDATION) {
        assert(lattice->sites[site_one].species == reaction.reactants[0]);
        assert(lattice->sites[site_two].species == SPECIES_EMPTY);

        relevant_react(lattice, update_function, site_one, std::optional<int> ());
    }


} //update_propensities() lattice

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::clear_site(Lattice *lattice, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props, 
                        int site, std::optional<int> ignore_neighbor, 
                        double &prop_sum, int &active_indices) {
    
    assert(site != SITE_GILLESPIE);
    // reset or initiate site combos
    for(uint32_t neigh = 0; neigh < lattice->numneigh[site]; neigh++ ) {
        int neighbor = lattice->idneigh[site][neigh];

        if(!ignore_neighbor || (ignore_neighbor && neighbor != ignore_neighbor)) {
            clear_site_helper(props, site, neighbor, prop_sum, active_indices);
        }
    }

    // self reactions 
    if(lattice->sites[site].species != SPECIES_EMPTY) {
        clear_site_helper(props, site, SITE_SELF_REACTION, prop_sum, active_indices);
    }
    

} // clear_site

/* ---------------------------------------------------------------------- */

// deal with active_indicies
void LatticeReactionNetwork::clear_site_helper(std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &props,
                        int site_one, int site_two, double &prop_sum,
                        int &active_indices) {

    std::string combo = make_string(site_one, site_two);
            
    // check if first time key has been added
    if(props.find(combo) == props.end()) {
        // key not found
        props[combo].reserve(reactions.size());
    }
    else {
        // already exists, clear vector to update
        prop_sum -= sum_row(combo, props);
        active_indices -= props[combo].size();

        props[combo].clear();
    }

} // clear_site_helper

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::relevant_react(Lattice *lattice, std::function<void(LatticeUpdate lattice_update)> update_function,
                                            int site, std::optional<int> ignore_neighbor) {

    // all reactions related to central site 
    std::vector<int> potential_reactions = dependents[lattice->sites[site].species]; 

    // compute and add new propensities 
    for(unsigned long int reaction_id = 0; reaction_id < static_cast<int> (potential_reactions.size()); reaction_id++ ) {

        LatticeReaction reaction = reactions[reaction_id]; 

        // single reactant (must be of type lattice)
        if(reaction.number_of_reactants == 1 && reaction.phase_reactants[0] == Phase::LATTICE) {

            // if reaction produces a solid product make sure on edge 
            if(reaction.type == Type::DESORPTION) {
                if(std::find(lattice->can_desorb.begin(), lattice->can_desorb.end(), site) != lattice->can_desorb.end()) {
                    // in can_desorb vector
                    
                    double new_propensity = compute_propensity(1, 0, reaction_id, lattice);

                    update_function(LatticeUpdate {
                    .index = reaction_id,
                    .propensity = new_propensity,
                    .site_one = site,
                    .site_two = SITE_SELF_REACTION}); 
                }
            }
            else {
                // oxidation / reduction / homogeneous solid with one reactant 
                double new_propensity = compute_propensity(1, 0, reaction_id, lattice, site);

                update_function(LatticeUpdate {
                    .index = reaction_id,
                    .propensity = new_propensity,
                    .site_one = site,
                    .site_two = SITE_SELF_REACTION}); 
            }

        } // single reactant
        // two reactants 
        else if(reaction.number_of_reactants == 2 && reaction.type != Type::ADSORPTION) {
            
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
                    
                    if(!ignore_neighbor || (ignore_neighbor && neighbor != ignore_neighbor)) {
                        if(lattice->sites[neighbor].species == reaction.reactants[other_reactant_id]) {

                            double new_propensity = compute_propensity(1, 1, reaction_id, lattice);

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

std::string LatticeReactionNetwork::make_string(int site_one, int site_two) {

    return (site_one < site_two) ? std::to_string(site_one) + "." + 
    std::to_string(site_two) : std::to_string(site_two) + "." + std::to_string(site_one);

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

TrajectoriesSql LatticeReactionNetwork::history_element_to_sql(
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
    
    SqlStatement<LatticeReactionSql> reaction_statement (reaction_network_database);
    SqlReader<LatticeReactionSql> reaction_reader (reaction_statement);

    // reaction_id is lifted so we can do a sanity check after the
    // loop.  Make sure size of reactions vector, last reaction_id and
    // metadata number_of_reactions are all the same

    // read in Reactions
    while(std::optional<LatticeReactionSql> maybe_reaction_row = reaction_reader.next()) {

        LatticeReactionSql reaction_row = maybe_reaction_row.value();

        LatticeReaction reaction;

        reaction.number_of_reactants = static_cast<uint8_t> (reaction_row.number_of_reactants);
        reaction.number_of_products = static_cast<uint8_t> (reaction_row.number_of_products);
        reaction.reactants[0] = reaction_row.reactant_1;
        reaction.reactants[1] = reaction_row.reactant_2;
        reaction.products [0]= reaction_row.product_1;
        reaction.products[1] = reaction_row.product_2;
        reaction.rate = reaction_row.rate;
        reaction.dG = reaction_row.dG;
        reaction.reorganization = reaction_row.reorganization;

        reaction.phase_reactants[0] = (reaction_row.phase_reactant_1 == 'L') ? Phase::LATTICE : Phase::SOLUTION;
        reaction.phase_reactants[1] = (reaction_row.phase_reactant_2 == 'L') ? Phase::LATTICE : Phase::SOLUTION;
        reaction.phase_products[0] = (reaction_row.phase_product_1 == 'L') ? Phase::LATTICE : Phase::SOLUTION;
        reaction.phase_products[1] = (reaction_row.phase_product_2 == 'L') ? Phase::LATTICE : Phase::SOLUTION;

        if(reaction_row.type == 'O') {
            reaction.type = Type::OXIDATION;
        }
        else if(reaction_row.type == 'R') {
            reaction.type = Type::REDUCTION;
        }
        else if (reaction_row.type == 'F') {
            reaction.type = Type::DIFFUSION;
        }
        else if (reaction_row.type == 'S') {
            reaction.type =  Type::HOMOGENEOUS_SOLID;
        }
        else if (reaction_row.type == 'E') {
            reaction.type =  Type::HOMOGENEOUS_ELYTE;
        }
        else if (reaction_row.type == 'D') {
            reaction.type =  Type::DESORPTION;
        }
        else if (reaction_row.type == 'A') {
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

        for ( int i = 0; i < reaction.number_of_reactants; i++ ) {
            int reactant_id = reaction.reactants[i];
            dependents[reactant_id].push_back(reaction_id);
        }
    }

};

/* ---------------------------------------------------------------------- */

void LatticeReactionNetwork::update_state(std::vector<int> &state, int reaction_index) {
    
    LatticeReaction reaction = reactions[reaction_index];

    for (int m = 0;
         m < reaction.number_of_products;
         m++) {
        if(reaction.phase_reactants[m] == Phase::SOLUTION) state[reactions[reaction_index].reactants[m]]--;
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
        int reactant_id = reaction.reactants[i];
        species_of_interest.push_back(reactant_id);
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

    if(reaction.type == Type::OXIDATION || reaction.type == Type::REDUCTION) {
         bool reduction_in = (reaction.type == Type::REDUCTION) ? true : false;
         get_marcus_rate_coefficient(reaction.dG, reaction.reorganization, g_e, lattice->get_maxz(), reduction_in);
    }
    

    double p;
    // one reactant
    if (reaction.number_of_reactants == 1)
        p = state[reaction.reactants[0]] * reaction.rate;


    // two reactants
    else {
        if (reaction.reactants[0] == reaction.reactants[1])
            p = factor_duplicate
                * factor_two
                * state[reaction.reactants[0]]
                * (state[reaction.reactants[0]] - 1)
                * reaction.rate;

        else
            p = factor_two
                * state[reaction.reactants[0]]
                * state[reaction.reactants[1]]
                * reaction.rate;
    }

    return p;
}

/* ---------------------------------------------------------------------- */

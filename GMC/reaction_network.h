#ifndef REACTION_NETWORK_H
#define REACTION_NETWORK_H

#include <stdint.h>
#include <vector>
#include <optional>
#include <mutex>
#include <memory>
#include <unordered_map>
#include <cmath>
#include "../core/sql.h"
#include "../GMC/sql_types.h"
#include "../core/simulation.h"
#include "../LGMC/lattice.h"

class Reaction {
    public:
    // we assume that each reaction has zero, one or two reactants
    uint8_t number_of_reactants;
    uint8_t number_of_products;

    int reactants[2];
    int products[2];

    double rate;
};

enum Phase {LATTICE, SOLUTION};
enum Type {ADSORPTION, DESORPTION, HOMOGENEOUS_ELYTE, HOMOGENEOUS_SOLID, DIFFUSION,CHARGE_TRANSFER};


const double TUNNEL_COEF = 1.2;        // Electron tunneling coefficient, in A^-1 
const double TEMPERATURE = 298.15;     // In Kelvin
const double KB = 8.6173e-5;           // In eV/K
const double PLANCK = 4.1357e-15;      // In eV s


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


class LatticeReaction : public Reaction {
    public:
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

// parameters passed to the ReactionNetwork constructor
// by the dispatcher which are model specific
struct ReactionNetworkParameters {
};


class ReactionNetwork {
    public:
    std::vector<std::shared_ptr<Reaction>> reactions; // list of reactions (either Reaction or LatticeReactions)
    std::vector<int> initial_state; // initial state for all the simulations
    std::vector<double> initial_propensities; // initial propensities for all the reactions
    double factor_zero; // rate modifer for reactions with zero reactants
    double factor_two; // rate modifier for reactions with two reactants
    double factor_duplicate; // rate modifier for reactions of form A + A -> ...
     LGMC_NS::Lattice *initial_lattice;

    // maps species to the reactions which involve that species
    std::vector<std::vector<int>> dependents;

    ReactionNetwork(){};         // defualt constructor for inheritence

    ReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        ReactionNetworkParameters parameters);

    virtual void compute_dependents();

    virtual double compute_propensity(std::vector<int> &state,
        int reaction_index);

    virtual void update_state(std::vector<int> &state,
        int reaction_index);

    virtual void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction);

    // convert a history element as found a simulation to history
    // to a SQL type.
    TrajectoriesSql history_element_to_sql(
        int seed,
        HistoryElement history_element);

    // overwritten by base class
    virtual void init(SqlConnection &reaction_network_database);
    virtual void fill_reactions(SqlConnection &reaction_network_database);

    // for model compatibilty 
    void update_state(std::unordered_map<std::string,                     
                    std::vector< std::pair<double, int> > > &props,
                    std::vector<int> &state, int next_reaction, 
                    std::optional<int> site_one, std::optional<int> site_two, double prop_sum){assert(false);};

    void update_propensities(std::vector<int> &state, std::function<void(Update update)> update_function, 
                                    std::function<void(LatticeUpdate lattice_update)> lattice_update_function, 
                                    int next_reaction, std::optional<int> site_one, std::optional<int> site_two) {assert(false);};
};

ReactionNetwork::ReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     ReactionNetworkParameters)

    {
    initial_lattice = nullptr;

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
    factor_zero = factors_row.factor_zero;
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

    init(reaction_network_database);

};

void ReactionNetwork::init(SqlConnection &reaction_network_database) {
    // loading reactions
    // vectors are default initialized to empty.
    // it is "cleaner" to resize the default vector than to
    // drop it and reinitialize a new vector.

    fill_reactions(reaction_network_database);

    // computing initial propensities
    for (unsigned long int i = 0; i < initial_propensities.size(); i++) {
        initial_propensities[i] = compute_propensity(initial_state, i);
    }

    std::cerr << time_stamp() << "computing dependency graph...\n";
    compute_dependents();
    std::cerr << time_stamp() << "finished computing dependency graph\n";
}

void ReactionNetwork::compute_dependents() {
    // initializing dependency graph

    dependents.resize(initial_state.size());


    for ( unsigned int reaction_id = 0; reaction_id <  reactions.size(); reaction_id++ ) {
        std::shared_ptr<Reaction> reaction = reactions[reaction_id];

        for ( int i = 0; i < reaction->number_of_reactants; i++ ) {
            int reactant_id = reaction->reactants[i];
            dependents[reactant_id].push_back(reaction_id);
        }
    }

};


double ReactionNetwork::compute_propensity(std::vector<int> &state,
    int reaction_index) {

    std::shared_ptr<Reaction> reaction = reactions[reaction_index];

    double p;
    // zero reactants
    if (reaction->number_of_reactants == 0)
        p = factor_zero * reaction->rate;

    // one reactant
    else if (reaction->number_of_reactants == 1)
        p = state[reaction->reactants[0]] * reaction->rate;


    // two reactants
    else {
        if (reaction->reactants[0] == reaction->reactants[1])
            p = factor_duplicate
                * factor_two
                * state[reaction->reactants[0]]
                * (state[reaction->reactants[0]] - 1)
                * reaction->rate;

        else
            p = factor_two
                * state[reaction->reactants[0]]
                * state[reaction->reactants[1]]
                * reaction->rate;
    }

    return p;

};

void ReactionNetwork::update_state(std::vector<int> &state,
    int reaction_index) {

    for (int m = 0;
         m < reactions[reaction_index]->number_of_products;
         m++) {
        state[reactions[reaction_index]->reactants[m]]--;
    }

    for (int m = 0;
         m < reactions[reaction_index]->number_of_products;
         m++) {
        state[reactions[reaction_index]->products[m]]++;
    }

}


void ReactionNetwork::update_propensities(
    std::function<void(Update update)> update_function,
    std::vector<int> &state,
    int next_reaction
    ) {

    std::shared_ptr<Reaction> reaction = reactions[next_reaction];

    std::vector<int> species_of_interest;
    species_of_interest.reserve(4);

    for ( int i = 0; i < reaction->number_of_reactants; i++ ) {
        int reactant_id = reaction->reactants[i];
        species_of_interest.push_back(reactant_id);
    }


    for ( int j = 0; j < reaction->number_of_products; j++ ) {
        int product_id = reaction->products[j];
        species_of_interest.push_back(product_id);
    }


    for ( int species_id : species_of_interest ) {
        for ( unsigned int reaction_index : dependents[species_id] ) {

            double new_propensity = compute_propensity(state,
                reaction_index);

            update_function(Update {
                    .index = reaction_index,
                    .propensity = new_propensity});
        }
    }
}

void ReactionNetwork::fill_reactions(SqlConnection &reaction_network_database) {

    SqlStatement<ReactionSql> reaction_statement (reaction_network_database);
    SqlReader<ReactionSql> reaction_reader (reaction_statement);


    // reaction_id is lifted so we can do a sanity check after the
    // loop.  Make sure size of reactions vector, last reaction_id and
    // metadata number_of_reactions are all the same

    // read in Reactions
    while(std::optional<ReactionSql> maybe_reaction_row = reaction_reader.next()) {

        ReactionSql reaction_row = maybe_reaction_row.value();
        uint8_t number_of_reactants = reaction_row.number_of_reactants;
        uint8_t number_of_products = reaction_row.number_of_products;

        std::shared_ptr<Reaction> reaction = std::make_shared<Reaction>();
        reaction->number_of_reactants = number_of_reactants;
        reaction->number_of_products = number_of_products;
        reaction->reactants[0] = reaction_row.reactant_1;
        reaction->reactants[1] = reaction_row.reactant_2;
        reaction->products [0]= reaction_row.product_1;
        reaction->products[1] = reaction_row.product_2;
        reaction->rate = reaction_row.rate;

        reactions.push_back(reaction);
    }

}

TrajectoriesSql ReactionNetwork::history_element_to_sql(
    int seed,
    HistoryElement history_element) {
    return TrajectoriesSql {
        .seed = seed,
        .step = history_element.step,
        .reaction_id = history_element.reaction_id,
        .time = history_element.time
    };
}


class LatticeReactionNetwork : public ReactionNetwork{
    public:
    void init(SqlConnection &reaction_network_database);
    void fill_reactions(SqlConnection &reaction_network_database);
    void compute_dependents();
    double compute_propensity(std::vector<int> &state, int reaction_index);
    void update_propensities(std::function<void(Update update)> update_function,
                             std::vector<int> &state, int next_reaction);
    void update_state(std::vector<int> &state, int reaction_index);

    LatticeReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        ReactionNetworkParameters parameters);

};

LatticeReactionNetwork::LatticeReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     ReactionNetworkParameters)

    {
    ReactionNetwork();
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
    factor_zero = factors_row.factor_zero;
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

    init(reaction_network_database);

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

        std::shared_ptr<LatticeReaction> reaction = std::make_shared<LatticeReaction>();

        reaction->number_of_reactants = static_cast<uint8_t> (reaction_row.number_of_reactants);
        reaction->number_of_products = static_cast<uint8_t> (reaction_row.number_of_products);
        reaction->reactants[0] = reaction_row.reactant_1;
        reaction->reactants[1] = reaction_row.reactant_2;
        reaction->products [0]= reaction_row.product_1;
        reaction->products[1] = reaction_row.product_2;
        reaction->rate = reaction_row.rate;
        reaction->dG = reaction_row.dG;
        reaction->reorganization = reaction_row.reorganization;

        reaction->phase_reactants[0] = (reaction_row.phase_reactant_1 == 'L') ? Phase::LATTICE : Phase::SOLUTION;
        reaction->phase_reactants[1] = (reaction_row.phase_reactant_2 == 'L') ? Phase::LATTICE : Phase::SOLUTION;
        reaction->phase_products[0] = (reaction_row.phase_product_1 == 'L') ? Phase::LATTICE : Phase::SOLUTION;
        reaction->phase_products[1] = (reaction_row.phase_product_2 == 'L') ? Phase::LATTICE : Phase::SOLUTION;

        if(reaction_row.type == 'C') {
            reaction->type = Type::CHARGE_TRANSFER;
        }
        else if (reaction_row.type == 'F') {
            reaction->type = Type::DIFFUSION;
        }
        else if (reaction_row.type == 'S') {
            reaction->type = Type::HOMOGENEOUS_SOLID;
        }
        else if (reaction_row.type == 'E') {
            reaction->type = Type::HOMOGENEOUS_ELYTE;
        }
        else if (reaction_row.type == 'D') {
            reaction->type = Type::DESORPTION;
        }
        else if (reaction_row.type == 'A') {
            reaction->type = Type::ADSORPTION;
        }

        reactions.push_back(reaction);
    }
}

void LatticeReactionNetwork::init(SqlConnection &reaction_network_database) {
    // loading reactions
    // vectors are default initialized to empty.
    // it is "cleaner" to resize the default vector than to
    // drop it and reinitialize a new vector.

    fill_reactions(reaction_network_database);

    // computing initial propensities
    for (unsigned long int i = 0; i < initial_propensities.size(); i++) {
        initial_propensities[i] = compute_propensity(initial_state, i);
    }

    std::cerr << time_stamp() << "computing dependency graph...\n";
    compute_dependents();
    std::cerr << time_stamp() << "finished computing dependency graph\n";
}

void LatticeReactionNetwork::compute_dependents() {
    // initializing dependency graph

    dependents.resize(initial_state.size());

    for ( unsigned int reaction_id = 0; reaction_id <  reactions.size(); reaction_id++ ) {
        LatticeReaction *reaction = static_cast<LatticeReaction*> (reactions[reaction_id].get());

        for ( int i = 0; i < reaction->number_of_reactants; i++ ) {
            int reactant_id = reaction->reactants[i];
            dependents[reactant_id].push_back(reaction_id);
        }
    }

};

void LatticeReactionNetwork::update_state(std::vector<int> &state, int reaction_index) {
    
    LatticeReaction *reaction = static_cast<LatticeReaction *> (reactions[reaction_index].get());

    for (int m = 0;
         m < reactions[reaction_index]->number_of_products;
         m++) {
        if(reaction->phase_reactants[m] == Phase::SOLUTION) state[reactions[reaction_index]->reactants[m]]--;
    }

    for (int m = 0;
         m < reactions[reaction_index]->number_of_products;
         m++) {
        if(reaction->phase_products[m] == Phase::SOLUTION) state[reactions[reaction_index]->products[m]]++;
    }
};


void LatticeReactionNetwork::update_propensities(std::function<void(Update update)> update_function,
                                                 std::vector<int> &state, int next_reaction) {

    LatticeReaction *reaction = static_cast<LatticeReaction *> (reactions[next_reaction].get());

    std::vector<int> species_of_interest;
    species_of_interest.reserve(4);

    for ( int i = 0; i < reaction->number_of_reactants; i++ ) {
        int reactant_id = reaction->reactants[i];
        species_of_interest.push_back(reactant_id);
    }


    for ( int j = 0; j < reaction->number_of_products; j++ ) {
        // make sure product is in the solution
        if(reaction->phase_products[j] == Phase::SOLUTION) {
            int product_id = reaction->products[j];
            species_of_interest.push_back(product_id);
        }
        
    }

    for ( int species_id : species_of_interest ) {
        for ( unsigned int reaction_index : dependents[species_id] ) {

            double new_propensity = compute_propensity(state, reaction_index);

            update_function(Update {
                    .index = reaction_index,
                    .propensity = new_propensity});
        }
    }
}

double LatticeReactionNetwork::compute_propensity(std::vector<int> &state, int reaction_index) {
    LatticeReaction *reaction = static_cast<LatticeReaction *> (reactions[reaction_index].get());

    double p;
    // zero reactants
    if (reaction->number_of_reactants == 0)
        p = factor_zero * reaction->rate;

    // one reactant
    else if (reaction->number_of_reactants == 1)
        p = state[reaction->reactants[0]] * reaction->rate;


    // two reactants
    else {
        if (reaction->reactants[0] == reaction->reactants[1])
            p = factor_duplicate
                * factor_two
                * state[reaction->reactants[0]]
                * (state[reaction->reactants[0]] - 1)
                * reaction->rate;

        else
            p = factor_two
                * state[reaction->reactants[0]]
                * state[reaction->reactants[1]]
                * reaction->rate;
    }

    return p;
}
#endif
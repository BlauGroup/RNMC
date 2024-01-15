#ifndef RNMC_GILLESPIE_REACTION_NETWORK_H
#define RNMC_GILLESPIE_REACTION_NETWORK_H

#include "reaction_network.h"

// parameters passed to the ReactionNetwork constructor
// by the dispatcher which are model specific
struct ReactionNetworkParameters {
    bool isCheckpoint;
};

struct GillespieReaction {
    // we assume that each reaction has zero, one or two reactants
    uint8_t number_of_reactants;
    uint8_t number_of_products;

    int reactants[2];
    int products[2];

    double rate;
};

class GillespieReactionNetwork : public ReactionNetwork<GillespieReaction> {
public:
    uint8_t energy_budget = 0;

    GillespieReactionNetwork(
        SqlConnection &reaction_network_database,
        SqlConnection &initial_state_database,
        ReactionNetworkParameters parameters);

    void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction);

    void checkpoint(SqlReader<ReactionNetworkReadStateSql> state_reader, 
        SqlReader<ReadCutoffSql> cutoff_reader, 
        SqlReader<ReactionNetworkReadTrajectoriesSql> trajectory_reader, 
        std::map<int, std::vector<int>> &temp_seed_state_map, 
        std::map<int, int> &temp_seed_step_map, 
        SeedQueue &temp_seed_queue, 
        std::map<int, double> &temp_seed_time_map, 
        ReactionNetwork &model);

    void store_checkpoint(std::vector<ReactionNetworkStateHistoryElement> 
        &state_packet, std::vector<int> &state,
        unsigned long int &seed, int step, double time, 
        std::vector<CutoffHistoryElement> &cutoff_packet);
        
    WriteCutoffSql cutoff_history_element_to_sql(
        int seed,
        CutoffHistoryElement cutoff_history_element);

};

GillespieReactionNetwork::GillespieReactionNetwork(
     SqlConnection &reaction_network_database,
     SqlConnection &initial_state_database,
     ReactionNetworkParameters parameters)
    {

    isCheckpoint = parameters.isCheckpoint;

    // collecting reaction network metadata
    SqlStatement<MetadataSql> metadata_statement (reaction_network_database);
    SqlReader<MetadataSql> metadata_reader (metadata_statement);

    std::optional<MetadataSql> maybe_metadata_row = metadata_reader.next();

    if (! maybe_metadata_row.has_value()) {
        std::cerr << time::time_stamp()
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

    // loading reactions
    // vectors are default initialized to empty.
    // it is "cleaner" to resize the default vector than to
    // drop it and reinitialize a new vector.
    reactions.resize(metadata_row.number_of_reactions);
    initial_propensities.resize(metadata_row.number_of_reactions);

    SqlStatement<ReactionSql> reaction_statement (reaction_network_database);
    SqlReader<ReactionSql> reaction_reader (reaction_statement);

    // reaction_id is lifted so we can do a sanity check after the
    // loop.  Make sure size of reactions vector, last reaction_id and
    // metadata number_of_reactions are all the same

    unsigned long int reaction_id = 0;

    while(std::optional<ReactionSql> maybe_reaction_row = reaction_reader.next()) {

        ReactionSql reaction_row = maybe_reaction_row.value();
        uint8_t number_of_reactants = reaction_row.number_of_reactants;
        uint8_t number_of_products = reaction_row.number_of_products;
        reaction_id = reaction_row.reaction_id;

        GillespieReaction reaction = {
            .number_of_reactants = number_of_reactants,
            .number_of_products = number_of_products,
            .reactants = { reaction_row.reactant_1, reaction_row.reactant_2 },
            .products = { reaction_row.product_1, reaction_row.product_2},
            .rate = reaction_row.rate
        };
        reactions[reaction_id] = reaction;
    }
     
    std::cerr << time::time_stamp() << "computing dependency graph...\n";
    compute_dependents();
    std::cerr << time::time_stamp() << "finished computing dependency graph\n";

} // GillespieReactionNetwork()

/*---------------------------------------------------------------------------*/

void GillespieReactionNetwork::update_propensities(
    std::function<void(Update update)> update_function,
    std::vector<int> &state,
    int next_reaction) {

    GillespieReaction &reaction = reactions[next_reaction];

    std::vector<int> species_of_interest;
    species_of_interest.reserve(4);

    for ( int i = 0; i < reaction.number_of_reactants; i++ ) {
        int reactant_id = reaction.reactants[i];
        species_of_interest.push_back(reactant_id);
    }

    for ( int j = 0; j < reaction.number_of_products; j++ ) {
        int product_id = reaction.products[j];
        species_of_interest.push_back(product_id);
    }

    for ( int species_id : species_of_interest ) {
        for ( unsigned int reaction_index : dependents[species_id] ) {

            double new_propensity = compute_propensity(
                state,
                reaction_index);

            update_function(Update {
                    .index = reaction_index,
                    .propensity = new_propensity});
        }
    }
} // update_propensities()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void GillespieReactionNetwork::checkpoint(SqlReader<ReactionNetworkReadStateSql> state_reader, 
                                    SqlReader<ReadCutoffSql> cutoff_reader, 
                                    SqlReader<ReactionNetworkReadTrajectoriesSql> trajectory_reader, 
                                    std::map<int, std::vector<int>> &temp_seed_state_map, 
                                    std::map<int, int> &temp_seed_step_map, 
                                    SeedQueue &temp_seed_queue, 
                                    std::map<int, double> &temp_seed_time_map, 
                                    ReactionNetwork<Reaction> &model) {
    
    bool read_interrupt_states = false;
    std::vector<int> default_state = model.initial_state;

    while (std::optional<unsigned long int> maybe_seed =
               temp_seed_queue.get_seed()){
                unsigned long int seed = maybe_seed.value();
                temp_seed_state_map.insert(std::make_pair(seed, default_state));
    }

    while (std::optional<ReadCutoffSql> maybe_cutoff_row = cutoff_reader.next()){
            ReadCutoffSql cutoff_row = maybe_cutoff_row.value();

            temp_seed_step_map[cutoff_row.seed] = cutoff_row.step;
            temp_seed_time_map[cutoff_row.seed] = cutoff_row.time;
    }

    while (std::optional<ReactionNetworkReadStateSql> maybe_state_row = state_reader.next()){
        read_interrupt_states = true;

        ReactionNetworkReadStateSql state_row = maybe_state_row.value();
        temp_seed_state_map[state_row.seed][state_row.species_id] = state_row.count;
    }

    if(!read_interrupt_states && isCheckpoint) {
        while (std::optional<ReactionNetworkReadTrajectoriesSql> maybe_trajectory_row = trajectory_reader.next()) {

            ReactionNetworkReadTrajectoriesSql trajectory_row = maybe_trajectory_row.value();
            
            Reaction reaction = model.reactions[trajectory_row.reaction_id];
            // update reactants
            for (int i = 0; i < reaction.number_of_reactants; i++) {
                temp_seed_state_map[trajectory_row.seed][reaction.reactants[i]] = 
                temp_seed_state_map[trajectory_row.seed][reaction.reactants[i]] - 1;
            }
            // update products
            for (int i = 0; i < reaction.number_of_products; i++) {
                temp_seed_state_map[trajectory_row.seed][reaction.products[i]] = 
                temp_seed_state_map[trajectory_row.seed][reaction.products[i]] + 1;
            }

            if (trajectory_row.step > temp_seed_step_map[trajectory_row.seed]) {
                temp_seed_step_map[trajectory_row.seed] = trajectory_row.step;
                temp_seed_time_map[trajectory_row.seed] = trajectory_row.time;
            }
        }
    }
}

/*---------------------------------------------------------------------------*/

template <typename Reaction>
void GillespieReactionNetwork::store_checkpoint(std::vector<ReactionNetworkStateHistoryElement> 
    &state_packet, std::vector<int> &state,
    unsigned long int &seed, int step, double time, 
    std::vector<CutoffHistoryElement> &cutoff_packet) {
    
    // state information
    for (unsigned int i = 0; i < state.size(); i++) {
        state_packet.push_back(ReactionNetworkStateHistoryElement{
            .seed = seed,
            .species_id = static_cast<int>(i),
            .count = state[i]
        });
    }

    // cutoff information
    cutoff_packet.push_back(CutoffHistoryElement {
        .seed = seed,
        .step = step,
        .time = time
    });
} // store_state_history()

/*---------------------------------------------------------------------------*/

template <typename Reaction>
WriteCutoffSql ReactionNetwork<Reaction>::cutoff_history_element_to_sql(
    int seed, CutoffHistoryElement cutoff_history_element) {
        return WriteCutoffSql {
            .seed = seed,
            .step = cutoff_history_element.step,
            .time = cutoff_history_element.time
        };
} // cutoff_history_element_to_sql()

#endif
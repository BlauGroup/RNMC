#include "nano_particle.h"

NanoParticle::NanoParticle() {
};

NanoParticle::NanoParticle(
    SqlConnection &nano_particle_database,
    SqlConnection &initial_state_database,
    NanoParticleParameters parameters) {

    isCheckpoint = parameters.isCheckpoint;

    // sql statements
    SqlStatement<SpeciesSql> species_statement(nano_particle_database);
    SqlStatement<SiteSql> site_statement(nano_particle_database);
    SqlStatement<InteractionSql> interactions_statement(nano_particle_database);
    SqlStatement<NanoMetadataSql> metadata_statement(nano_particle_database);
    SqlStatement<NanoFactorsSql> factors_statement(initial_state_database);
    SqlStatement<NanoInitialStateSql> initial_state_statement(initial_state_database);

    // sql readers
    SqlReader<SpeciesSql> species_reader(species_statement);
    SqlReader<SiteSql> site_reader(site_statement);
    SqlReader<InteractionSql> interactions_reader(interactions_statement);
    SqlReader<NanoMetadataSql> metadata_reader(metadata_statement);
    SqlReader<NanoFactorsSql> factors_reader(factors_statement);
    SqlReader<NanoInitialStateSql> initial_state_reader(initial_state_statement);

    // extracting metadata
    std::optional<NanoMetadataSql> maybe_metadata_row =
        metadata_reader.next();

    if (! maybe_metadata_row.has_value()) {
        std::cerr << time::time_stamp()
                  << "no metadata row\n";

        std::abort();
    }

    NanoMetadataSql metadata_row = maybe_metadata_row.value();

    // extracting factors
    std::optional<NanoFactorsSql> maybe_factor_row =
        factors_reader.next();

    if (! maybe_factor_row.has_value()) {
        std::cerr << time::time_stamp()
                  << "no factor row\n";

        std::abort();
    }

    NanoFactorsSql factor_row = maybe_factor_row.value();

    one_site_interaction_factor = factor_row.one_site_interaction_factor;
    two_site_interaction_factor = factor_row.two_site_interaction_factor;
    interaction_radius_bound = factor_row.interaction_radius_bound;

    if ( factor_row.distance_factor_type == "linear" ) {
        distance_factor_function = [=](double distance) {
            return 1 - ( distance / interaction_radius_bound ); };

    } else if ( factor_row.distance_factor_type == "inverse_cubic" ) {
        distance_factor_function = [](double distance) {
            return  1 / ( pow(distance,6)); };

    } else {
        std::cerr << time::time_stamp()
                  << "unexpected distance_factor_type: "
                  << factor_row.distance_factor_type << '\n'
                  << "expecting linear or inverse_cubic" << '\n';

       std::abort();
    }

    // initializing degrees of freedom
    degrees_of_freedom.resize(metadata_row.number_of_species);
    while(std::optional<SpeciesSql> maybe_species_row =
          species_reader.next()) {
        SpeciesSql species_row = maybe_species_row.value();

        degrees_of_freedom[species_row.species_id] =
            species_row.degrees_of_freedom;
    }

    // initializing sites
    sites.resize(metadata_row.number_of_sites);
    site_reaction_dependency.resize(metadata_row.number_of_sites);
    // for ( unsigned int i = 0; i < site_reaction_dependency.size(); i++) {
    //     site_reaction_dependency[i].resize(0);
    // }

    while(std::optional<SiteSql> maybe_site_row =
          site_reader.next()) {

        SiteSql site_row = maybe_site_row.value();
        sites[site_row.site_id] = {
            .x = site_row.x,
            .y = site_row.y,
            .z = site_row.z,
            .species_id = (int) site_row.species_id };
    }

    // initialize interactions
    // interactions.resize(metadata_row.number_of_interactions);
    int num_species = metadata_row.number_of_species;
    int interaction_counter = 0;
    int num_states = 0; // Keep track of number of states, so axis 2 and 3 in interaction_map can be resized
    while(std::optional<InteractionSql> maybe_interaction_row =
          interactions_reader.next()) {

        InteractionSql interaction_row = maybe_interaction_row.value();

        Interaction interaction = Interaction {
                .interaction_id  = interaction_counter,
                .number_of_sites = interaction_row.number_of_sites,
                .species_id      = { interaction_row.species_id_1, interaction_row.species_id_2},
                .left_state      = { interaction_row.left_state_1, interaction_row.left_state_2},
                .right_state     = { interaction_row.right_state_1, interaction_row.right_state_2},
                .rate            = interaction_row.rate
            };

        all_interactions.push_back(interaction);
        if (interaction_row.number_of_sites == 1) {
            one_site_interactions.push_back(interaction);
        } else if (interaction_row.number_of_sites == 2) {
            two_site_interactions.push_back(interaction);
        }

        if (num_states < interaction_row.left_state_1) {
            // Keep track of the max number of states
            num_states = interaction_row.left_state_1;
        }

        // Increment the interaction counter
        interaction_counter++;
    }
    // Increment the state counter, since this value will be off by 1
    num_states++;

    // Resize interaction_maps
    one_site_interactions_map.resize(num_species);
    for (int i = 0; i < num_species; i++) {
        one_site_interactions_map[i].resize(num_states);
    }

    two_site_interactions_map.resize(num_species);
    for (int i = 0; i < num_species; i++) {
        two_site_interactions_map[i].resize(num_species);
        for (int j = 0; j < num_species; j++) {
            two_site_interactions_map[i][j].resize(num_states);
            for (int k = 0; k < num_states; k++) {
                two_site_interactions_map[i][j][k].resize(num_states);
            }
        }
    }

    // Put interactions into interaction_map
    for (unsigned int i = 0; i < one_site_interactions.size(); i++) {
        Interaction interaction = one_site_interactions[i];
        one_site_interactions_map[interaction.species_id[0]][interaction.left_state[0]].push_back(interaction);
    }

    for (unsigned int i = 0; i < two_site_interactions.size(); i++) {
        Interaction interaction = two_site_interactions[i];
        two_site_interactions_map[interaction.species_id[0]][interaction.species_id[1]][interaction.left_state[0]][interaction.left_state[1]].push_back(interaction);
    }

    // initialize initial_state
    initial_state.resize(metadata_row.number_of_sites);

    while(std::optional<NanoInitialStateSql> maybe_initial_state_row =
          initial_state_reader.next()) {
        NanoInitialStateSql initial_state_row = maybe_initial_state_row.value();
        initial_state[initial_state_row.site_id] = initial_state_row.degree_of_freedom;
    }

    // Pre-compute the distance matrix so that it doesn't need to be computed multiple times
    compute_distance_matrix();
} // NanoParticle()

/* ---------------------------------------------------------------------- */

double NanoParticle::site_distance_squared(NanoSite s1, NanoSite s2) {
    double x_diff = s1.x - s2.x;
    double y_diff = s1.y - s2.y;
    double z_diff = s1.z - s2.z;

    return ( x_diff * x_diff +
             y_diff * y_diff +
             z_diff * z_diff );

} // site_distance_squared()

/* ---------------------------------------------------------------------- */

void NanoParticle::compute_reactions(
    const std::vector<int> &state,
    std::vector<NanoReaction> &reactions,    
    std::vector<std::set<int>> &site_reaction_dependency) {
    
    int reaction_count = 0;
    for (unsigned int site_id_0 = 0; site_id_0 < sites.size(); site_id_0++) {
        // Add one site interactions
        int site_0_state = state[site_id_0];
        int site_0_species_id = sites[site_id_0].species_id;
        std::vector<Interaction>* available_interactions = &one_site_interactions_map[site_0_species_id][site_0_state];
        for (unsigned int i = 0; i < available_interactions->size(); i++) {
            Interaction interaction = (*available_interactions)[i];
            NanoReaction reaction = NanoReaction {
                                .site_id = { (int) site_id_0, -1},
                                .interaction = interaction,
                                .rate = interaction.rate * one_site_interaction_factor
                                };
            reactions.push_back(reaction);
            site_reaction_dependency[site_id_0].insert(reaction_count);
            reaction_count++;
        }

        // Add two site interactions
        for (unsigned int site_id_1 = 0; site_id_1 < sites.size(); site_id_1++) {
            if ((unsigned) site_id_0 != site_id_1) {
                int site_1_state = state[site_id_1];
                int site_1_species_id = sites[site_id_1].species_id;
                double distance = distance_matrix[site_id_0][site_id_1];
                if (distance < interaction_radius_bound) {
                    // Add reactions where site 0 is the donor
                    std::vector<Interaction>* available_interactions = &two_site_interactions_map[site_0_species_id][site_1_species_id][site_0_state][site_1_state];
                    for (unsigned int i = 0; i < available_interactions->size(); i++){
                        Interaction interaction = (*available_interactions)[i];
                        NanoReaction reaction = NanoReaction {
                                                .site_id = { (int) site_id_0, (int) site_id_1 },
                                                .interaction = interaction,
                                                .rate = distance_factor_function(distance) * interaction.rate * two_site_interaction_factor
                                                };
                        reactions.push_back(reaction);
                        site_reaction_dependency[site_id_0].insert(reaction_count);
                        site_reaction_dependency[site_id_1].insert(reaction_count);
                        reaction_count++;
                    }

                    // Add reactions where site 1 is the donor
                    available_interactions = &two_site_interactions_map[site_1_species_id][site_0_species_id][site_1_state][site_0_state];
                    for (unsigned int i = 0; i < available_interactions->size(); i++){
                        Interaction interaction = (*available_interactions)[i];
                        NanoReaction reaction = NanoReaction {
                                                .site_id = { (int) site_id_1, (int) site_id_0 },
                                                .interaction = interaction,
                                                .rate = distance_factor_function(distance) * interaction.rate * two_site_interaction_factor
                                                };
                        reactions.push_back(reaction);
                        site_reaction_dependency[site_id_0].insert(reaction_count);
                        site_reaction_dependency[site_id_1].insert(reaction_count);
                        reaction_count++;
                    }
                }
            }
        }
    }
} // compute_reactions()

/* ---------------------------------------------------------------------- */

void NanoParticle::compute_distance_matrix() {
    distance_matrix.resize(sites.size());
    for (unsigned int i = 0; i < sites.size(); i++) {
        distance_matrix[i].resize(sites.size());
        for (unsigned int j = 0; j < sites.size(); j++) {
            distance_matrix[i][j] = std::sqrt(site_distance_squared(sites[i], sites[j]));
        }
    }
} // compute_distance_matrix()

/* ---------------------------------------------------------------------- */

void NanoParticle::update_state(
    std::vector<int> &state,
    NanoReaction reaction) {

    Interaction interaction = reaction.interaction;

    for (int k = 0; k < interaction.number_of_sites; k++) {
        if (state[reaction.site_id[k]] != interaction.left_state[k]) {
            // Ensure that the update is valid
            std::cerr << "State mismatch for site_id "
                      << reaction.site_id[k]
                      << "Expected state"
                      << interaction.left_state[k]
                      << ", found state "
                      << state[reaction.site_id[k]]
                      << "\n";
            raise(SIGINT);
        }
        state[reaction.site_id[k]] = interaction.right_state[k];
    }
} // update_state()

/* ---------------------------------------------------------------------- */

void NanoParticle::compute_new_reactions(
    const int site_0_id,
    const int other_site_id,
    const int site_0_state,
    const std::vector<int> &state,
    std::vector<NanoReaction> &new_reactions){
        
    int site_0_species_id = sites[site_0_id].species_id;

    // Add one site interactions
    std::vector<Interaction>* available_interactions = &one_site_interactions_map[site_0_species_id][site_0_state];
    for (unsigned int i = 0; i < available_interactions->size(); i++) {
        Interaction interaction = (*available_interactions)[i];
        NanoReaction new_reaction = NanoReaction {
                                .site_id = { (int) site_0_id, -1},
                                .interaction = interaction,
                                .rate = interaction.rate * one_site_interaction_factor
                                };
        new_reactions.push_back(new_reaction);
    }

    // Add two site interactions
    for (unsigned int site_1_id = 0; site_1_id < sites.size(); site_1_id++) {
        if ((unsigned) site_0_id != site_1_id) {
            int site_1_state = state[site_1_id];
            int site_1_species_id = sites[site_1_id].species_id;

            const double* distance = &distance_matrix[site_0_id][site_1_id];
            if (*distance < interaction_radius_bound) {
                // Add reactions where site 0 is the donor
                std::vector<Interaction>* available_interactions = &two_site_interactions_map[site_0_species_id][site_1_species_id][site_0_state][site_1_state];
                for (unsigned int i = 0; i < available_interactions -> size(); i++){
                    Interaction interaction = (*available_interactions)[i];
                    NanoReaction new_reaction = NanoReaction {
                                            .site_id = { (int) site_0_id, (int) site_1_id },
                                            .interaction = interaction,
                                            .rate = distance_factor_function(*distance) * interaction.rate * two_site_interaction_factor
                                            };
                    new_reactions.push_back(new_reaction);
                }

                // This if check is necessary so we don't doubly add reactions.
                // i.e. if our reaction which fired involves sites 11 and 22, we want to only add 11->22 and 22->11 once.
                // If this check isn't here, we add 11->22 and 22->11 twice
                if (site_1_id != (unsigned) other_site_id) {
                    // Add reactions where site 1 is the donor
                    available_interactions = &two_site_interactions_map[site_1_species_id][site_0_species_id][site_1_state][site_0_state];
                    for (unsigned int i = 0; i < available_interactions->size(); i++){
                        Interaction interaction = (*available_interactions)[i];
                        NanoReaction new_reaction = NanoReaction {
                                                .site_id = { (int) site_1_id, (int) site_0_id },
                                                .interaction = interaction,
                                                .rate = distance_factor_function(*distance) * interaction.rate * two_site_interaction_factor
                                                };
                        new_reactions.push_back(new_reaction);
                    }
                }
            }
        }
    }
} // compute_new_reactions()

/* ---------------------------------------------------------------------- */

void NanoParticle::update_reactions(
    const std::vector<int>& state,
    NanoReaction reaction,
    std::vector<std::set<int>>& current_site_reaction_dependency,
    std::vector<NanoReaction>& current_reactions) {

    // Compute the new reactions based on the new states
    std::vector<NanoReaction> new_reactions;
    const int site_0_id = reaction.site_id[0];
    const int site_0_state = state[site_0_id];
    const int site_1_id = reaction.site_id[1];
    const int site_1_state = state[site_1_id];
    compute_new_reactions(site_0_id, site_1_id, site_0_state, std::cref(state), std::ref(new_reactions));
    if (reaction.interaction.number_of_sites == 2) {
        compute_new_reactions(site_1_id, site_0_id, site_1_state, std::cref(state), std::ref(new_reactions));
    }

    // Find indexes of reactions to be removed
    std::set<int> reactions_to_remove;
    std::set<int>* reaction_dependency;
    for (int k = 0; k < reaction.interaction.number_of_sites; k++) {
        int site_id = reaction.site_id[k];
        reaction_dependency = &current_site_reaction_dependency[site_id];
        
        std::vector<std::vector<int>> dependencies_to_process;
        for (std::set<int>::iterator site_reaction_dependency_itr = reaction_dependency->begin();
            site_reaction_dependency_itr != reaction_dependency->end();
            site_reaction_dependency_itr++) {
            int reaction_id_to_remove = *site_reaction_dependency_itr;
            // Append this reaction index to a list of reactions to remove
            reactions_to_remove.insert(reaction_id_to_remove);

            // Need to remove this interaction from the second site if this is a two_site interaction
            NanoReaction* reaction_to_remove = &current_reactions[reaction_id_to_remove];
            dependencies_to_process.push_back({reaction_to_remove->site_id[0], reaction_id_to_remove});
            if ((*reaction_to_remove).interaction.number_of_sites == 2) {
                dependencies_to_process.push_back({reaction_to_remove->site_id[1], reaction_id_to_remove});
            }
        }

        for (unsigned int i = 0; i < dependencies_to_process.size(); i++){
            std::vector<int>* dependency_to_process = &dependencies_to_process[i];
            current_site_reaction_dependency[(*dependency_to_process)[0]].erase((*dependency_to_process)[1]);
        }
    }
    
    // Add the new reactions to the current_reactions vector
    int n_reactions_to_remove = reactions_to_remove.size();
    std::set<int>::iterator reactions_to_remove_itr = reactions_to_remove.begin();
    for (unsigned int i = 0; i < new_reactions.size(); i++) {
        NanoReaction* new_reaction = &new_reactions[i];
        if (i < (unsigned) n_reactions_to_remove) {
            // Assign the new reaction to a index belonging to a reaction to remove.
            // Since the reaction is going to be deleted anyways, this is safe.
            // Additionally, this avoids additional copy operations
            current_reactions[*reactions_to_remove_itr] = *new_reaction;
            for (int k = 0; k < (*new_reaction).interaction.number_of_sites; k++) {
                current_site_reaction_dependency[new_reaction->site_id[k]].insert(*reactions_to_remove_itr);
            }
            // reactions_to_remove.erase(reactions_to_remove_itr++);
            reactions_to_remove_itr++;
        } else {
            // If the number of new reactions to be added is larger than the number of reactions to remove,
            // just append the excess reactions to the end of the current_reactions vector
            current_reactions.push_back(*new_reaction);
            for (int k = 0; k < (*new_reaction).interaction.number_of_sites; k++) {
                current_site_reaction_dependency[new_reaction->site_id[k]].insert(current_reactions.size()-1);
            }
        }
    }

    int net_change_in_num_reactions = new_reactions.size() - reactions_to_remove.size();
    if (reactions_to_remove_itr != reactions_to_remove.end() && net_change_in_num_reactions < 0){
        // Don't need to do anything with the reactions_to_remove_itr, since it already points to the
        // last reaction that has yet to be removed (a result of the last code block)

        // Create a new variable keeping track of which reaction (within current_reactions) we are checking.
        // This will start at the end of current_reactions.
        int reaction_idx_to_move = (int) current_reactions.size() - 1;
        int reactions_moved = 0;
        while (reactions_to_remove_itr != reactions_to_remove.end() &&
               reactions_moved > net_change_in_num_reactions &&
               reaction_idx_to_move >= *reactions_to_remove_itr) {
            // Check if the reaction to be moved is in the list of reactions to remove
            auto result = reactions_to_remove.find(reaction_idx_to_move);
            if (result != reactions_to_remove.end()) {
                // Reaction to be moved is going to be removed, no point in removing it
                reaction_idx_to_move--;
            } else {
                NanoReaction reaction_to_move = current_reactions[reaction_idx_to_move];
                current_reactions[*reactions_to_remove_itr] = reaction_to_move;

                // Find the reaction that was moved in the site reaction dependency vector and remap it
                for (int k = 0; k < reaction_to_move.interaction.number_of_sites; k++) {
                    int site_id = reaction_to_move.site_id[k];
                    // This number should always be in the list, else something has gone wrong
                    int num_reactions_deleted = current_site_reaction_dependency[site_id].erase(reaction_idx_to_move);

                    if (num_reactions_deleted == 1) {
                        current_site_reaction_dependency[site_id].insert(*reactions_to_remove_itr);
                    } else if (num_reactions_deleted == 0) {
                        std::cerr << "Could not find reaction "
                                    << reaction_idx_to_move
                                    << " in the site reaction dependency map for site "
                                    << site_id
                                    << "\n";
                        raise(SIGINT);
                    } else {
                        std::cerr << "Encountered multiple occurances of "
                                    << reaction_idx_to_move
                                    << " in the site reaction dependency map for site "
                                    << site_id
                                    << "\n";
                        raise(SIGINT);
                    }
                }
                reaction_idx_to_move--;
                reactions_to_remove_itr++;
                reactions_moved++;
            }
        }
        current_reactions.resize(current_reactions.size() + net_change_in_num_reactions);
    }

} // update_reactions()

/* ---------------------------------------------------------------------- */

NanoWriteTrajectoriesSql NanoParticle::history_element_to_sql(
    int seed,
    NanoTrajectoryHistoryElement history_element) {

    NanoReaction reaction = history_element.reaction;
    return NanoWriteTrajectoriesSql {
        .seed = seed,
        .step = history_element.step,
        .time = history_element.time,
        .site_id_1 = reaction.site_id[0],
        .site_id_2 = reaction.site_id[1],
        .interaction_id = reaction.interaction.interaction_id
    };
} // history_element_to_sql()

/* ---------------------------------------------------------------------- */

NanoWriteStateSql NanoParticle::state_history_element_to_sql(
    int seed,
    NanoStateHistoryElement state_history_element) {

    return NanoWriteStateSql {
        .seed = seed,
        .site_id = state_history_element.site_id,
        .degree_of_freedom = state_history_element.degree_of_freedom
    };
} // state_history_element_to_sql()

/* ---------------------------------------------------------------------- */

WriteCutoffSql NanoParticle::cutoff_history_element_to_sql(
    int seed,
    CutoffHistoryElement cutoff_history_element) {
        return WriteCutoffSql {
            .seed = seed,
            .step = cutoff_history_element.step,
            .time = cutoff_history_element.time
        };
} // cutoff_history_element_to_sql()

/* ---------------------------------------------------------------------- */

void NanoParticle::checkpoint(SqlReader<NanoReadStateSql> state_reader, 
                    SqlReader<ReadCutoffSql> cutoff_reader, 
                    SqlReader<NanoReadTrajectoriesSql> trajectory_reader, 
                    std::map<int, std::vector<int>> &temp_seed_state_map, 
                    std::map<int, int> &temp_seed_step_map, 
                    SeedQueue &temp_seed_queue, 
                    std::map<int, double> &temp_seed_time_map, 
                    NanoParticle &model) {
    
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

    // try reading from state
    while (std::optional<NanoReadStateSql> maybe_state_row = state_reader.next()){
        read_interrupt_states = true;

        NanoReadStateSql state_row = maybe_state_row.value();
        temp_seed_state_map[state_row.seed][state_row.site_id] = state_row.degree_of_freedom;
    }

    // try reading from trajectory
    if(!read_interrupt_states) {
        while (std::optional<NanoReadTrajectoriesSql> maybe_trajectory_row = trajectory_reader.next()) {
                            // read_trajectory_states = true;

            NanoReadTrajectoriesSql trajectory_row = maybe_trajectory_row.value();
            
            Interaction *interaction = &model.all_interactions[trajectory_row.interaction_id];
            temp_seed_state_map[trajectory_row.seed][trajectory_row.site_id_1] = interaction->right_state[0];
            if (interaction->number_of_sites == 2) {
                temp_seed_state_map[trajectory_row.seed][trajectory_row.site_id_2] = interaction->right_state[1];
            }

            if (trajectory_row.step > temp_seed_step_map[trajectory_row.seed]) {
                temp_seed_step_map[trajectory_row.seed] = trajectory_row.step;
                temp_seed_time_map[trajectory_row.seed] = trajectory_row.time;
            }
        }
    }
} // checkpoint()

/* ---------------------------------------------------------------------- */

void NanoParticle::store_checkpoint(std::vector<NanoStateHistoryElement> &state_packet,
    std::vector<int> &state,
    unsigned long int &seed, int step, double time, 
    std::vector<CutoffHistoryElement> &cutoff_packet) {   

    // state information
    for (unsigned int i = 0; i < state.size(); i++) {
                state_packet.push_back(NanoStateHistoryElement{
                    .seed = seed,
                    .site_id = static_cast<int>(i),
                    .degree_of_freedom = state[i]
                });
    }

    // cutoff information
    cutoff_packet.push_back(CutoffHistoryElement {
        .seed = seed,
        .step = step,
        .time = time
    });
} // store_state_history()

/* ---------------------------------------------------------------------- */
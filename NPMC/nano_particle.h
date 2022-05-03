#pragma once
#include "../core/sql.h"
#include "../core/simulation.h"
#include "sql_types.h"
#include <vector>
#include <cmath>
#include <functional>
#include <csignal>
#include <set>

struct Site {
    double x;
    double y;
    double z;
    int species_id;
};

double site_distance_squared(Site s1, Site s2) {
    double x_diff = s1.x - s2.x;
    double y_diff = s1.y - s2.y;
    double z_diff = s1.z - s2.z;

    return ( x_diff * x_diff +
             y_diff * y_diff +
             z_diff * z_diff );
}

// struct Interaction {
//     // either 1 site or two site interaction
//     int interaction_id;
//     int number_of_sites;
//     int species_id[2];
//     int left_state[2];
//     int right_state[2];
//
//     // the units of rate depend on number_of_sites. If number_of_sites = 1, then
//     // rate has units 1 / s. If number of sites = 2, then rate has units 1 / s m^6.
//     double rate;
// };

// For this nano particle simulator, a reaction consists
// of upto two sites and the interaction id.
// if it is an internal interaction, site_id_2 will be -1
// In a reaction, the sites must be within the interaction radius bound.


// struct Reaction {
//     int site_id[2];
//     Interaction interaction;
//
//     // rate has units 1 / s
//     double rate;
// };

struct NanoParticleParameters {};

struct NanoParticle {
    // maps a species index to the number of degrees of freedom
    std::vector<int> degrees_of_freedom;

    // maps site index to site data
    std::vector<Site> sites;

    // 2D vector representing the pairwise distance between two sites
    std::vector<std::vector<double>> distance_matrix;

    // maps site ids to sites which are within a cutoff distance
    std::vector<std::vector<int>> site_dependency;

    // maps site ids to sit
    std::vector<Reaction> initial_reactions;

    // Remove this since we don't have concrete reaction_ids
    // maps site ids to reaction ids involving the site.
    std::vector<std::set<int>> site_reaction_dependency;

    // maps interaction index to interaction data
    std::vector<Interaction> one_site_interactions;
    std::vector<Interaction> two_site_interactions;

    // maps site_0_index, site_0_state to interaction data
    std::vector<std::vector<std::vector<Interaction>>> one_site_interactions_map;
    // maps site_0_index, site_1_index, site_0_state, site_1_state to interaction data
    std::vector<std::vector<std::vector<std::vector<std::vector<Interaction>>>>> two_site_interactions_map;

    // initial state of the simulations.
    // initial_state[i] is a local degree of freedom
    // from the species at site i.
    std::vector<int> initial_state;
    std::vector<int> current_state;

    // std::vector<double> initial_propensities;

    // list mapping reaction_ids to reactions
    std::vector<Reaction> reactions;

    double one_site_interaction_factor;
    double two_site_interaction_factor;
    double interaction_radius_bound;

    std::function<double(double)> distance_factor_function;

    // constructor
    NanoParticle(
        SqlConnection &nano_particle_database,
        SqlConnection &initial_state_database,
        NanoParticleParameters
        );

    // maps a site index to the indices of its neighbors
    // within the spatial decay radius
    std::vector<std::vector<int>> compute_site_neighbors();

    void compute_reactions();
    void compute_distance_matrix();
    void find_dependency();

    double compute_propensity(
        std::vector<int> &state,
        Reaction reaction);

    void update_state(
        std::vector<int> &state,
        Reaction reaction);

    // updates are passed directly to the solver, but the model
    // classes don't have a reference to the solver (otherwise,
    // the would need to be parameterized over a Solver type, which
    // feels akward). The way around this is that the model class computes
    // all the propensity updates required, and then passes them into an
    // update_function lambda which is generated by the Simulation class from
    // core.

    void update_reactions(
        std::vector<int> &state,
        std::vector<std::set<int>> &current_site_reaction_dependency,
        std::vector<Reaction> &current_reactions,
        Reaction reaction);

    // void update_propensities(
    //     std::function<void(Update update)> update_function,
    //     std::vector<int> &state,
    //     Reaction reaction
    //     );

    // convert a history element as found a simulation to history
    // to a SQL type.
    TrajectoriesSql history_element_to_sql(
        int seed,
        HistoryElement history_element);

};

NanoParticle::NanoParticle(
    SqlConnection &nano_particle_database,
    SqlConnection &initial_state_database,
    NanoParticleParameters
    ) {

    // sql statements
    SqlStatement<SpeciesSql> species_statement(nano_particle_database);
    SqlStatement<SiteSql> site_statement(nano_particle_database);
    SqlStatement<InteractionSql> interactions_statement(nano_particle_database);
    SqlStatement<MetadataSql> metadata_statement(nano_particle_database);
    SqlStatement<FactorsSql> factors_statement(initial_state_database);
    SqlStatement<InitialStateSql> initial_state_statement(initial_state_database);

    // sql readers
    SqlReader<SpeciesSql> species_reader(species_statement);
    SqlReader<SiteSql> site_reader(site_statement);
    SqlReader<InteractionSql> interactions_reader(interactions_statement);
    SqlReader<MetadataSql> metadata_reader(metadata_statement);
    SqlReader<FactorsSql> factors_reader(factors_statement);
    SqlReader<InitialStateSql> initial_state_reader(initial_state_statement);

    // extracting metadata
    std::optional<MetadataSql> maybe_metadata_row =
        metadata_reader.next();

    if (! maybe_metadata_row.has_value()) {
        std::cerr << time_stamp()
                  << "no metadata row\n";

        std::abort();
    }

    MetadataSql metadata_row = maybe_metadata_row.value();


    // extracting factors
    std::optional<FactorsSql> maybe_factor_row =
        factors_reader.next();

    if (! maybe_factor_row.has_value()) {
        std::cerr << time_stamp()
                  << "no factor row\n";

        std::abort();
    }

    FactorsSql factor_row = maybe_factor_row.value();

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
        std::cerr << time_stamp()
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

        if (interaction_row.number_of_sites == 1) {
            one_site_interactions.push_back( Interaction {
                .interaction_id  = interaction_counter,
                .number_of_sites = interaction_row.number_of_sites,
                .species_id      = { interaction_row.species_id_1, interaction_row.species_id_2},
                .left_state      = { interaction_row.left_state_1, interaction_row.left_state_2},
                .right_state     = { interaction_row.right_state_1, interaction_row.right_state_2},
                .rate            = interaction_row.rate
            });
        } else if (interaction_row.number_of_sites == 2) {
          two_site_interactions.push_back( Interaction {
              .interaction_id  = interaction_counter,
              .number_of_sites = interaction_row.number_of_sites,
              .species_id      = { interaction_row.species_id_1, interaction_row.species_id_2},
              .left_state      = { interaction_row.left_state_1, interaction_row.left_state_2},
              .right_state     = { interaction_row.right_state_1, interaction_row.right_state_2},
              .rate            = interaction_row.rate
          });
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

    while(std::optional<InitialStateSql> maybe_initial_state_row =
          initial_state_reader.next()) {
        InitialStateSql initial_state_row = maybe_initial_state_row.value();
        initial_state[initial_state_row.site_id] = initial_state_row.degree_of_freedom;
    }

    // Here we want to avoid computing reactions, so we only create a dependency list.
    // Also pre-compute the distance matrix so that it doesn't need to be computed multiple times
    compute_distance_matrix();
    // site_dependency = compute_site_neighbors();

    // Setup current_state by copying from initial_state
    current_state = initial_state;
    int reaction_count = 0;

    // raise(SIGINT);
    for (unsigned int site_id_0 = 0; site_id_0 < sites.size(); site_id_0++) {
        // Add one site interactions
        int site_0_state = current_state[site_id_0];
        int site_0_species_id = sites[site_id_0].species_id;
        std::vector<Interaction> available_interactions = one_site_interactions_map[site_0_species_id][site_0_state];
        for (unsigned int i = 0; i < available_interactions.size(); i++) {
            Interaction interaction = available_interactions[i];
            Reaction reaction = Reaction {
                                .site_id = { (int) site_id_0, -1},
                                .interaction = interaction,
                                .rate = interaction.rate * one_site_interaction_factor
                                };
            initial_reactions.push_back(reaction);
            site_reaction_dependency[site_id_0].insert(reaction_count);
            reaction_count++;
        }

        // Add two site interactions
        for (unsigned int site_id_1 = 0; site_id_1 < sites.size(); site_id_1++) {
            if (site_id_0 != (int) site_id_1) {
                int site_1_state = current_state[site_id_1];
                int site_1_species_id = sites[site_id_1].species_id;
                std::vector<Interaction> available_interactions = two_site_interactions_map[site_0_species_id][site_1_species_id][site_0_state][site_1_state];
                double distance = distance_matrix[site_id_0][site_id_1];
                if (distance <= interaction_radius_bound) {
                    for (unsigned int i = 0; i < available_interactions.size(); i++){
                        Interaction interaction = available_interactions[i];
                        Reaction reaction = Reaction {
                                            .site_id = { (int) site_id_0, (int) site_id_1 },
                                            .interaction = interaction,
                                            .rate = distance_factor_function(distance) * interaction.rate * two_site_interaction_factor
                                            };
                        initial_reactions.push_back(reaction);
                        site_reaction_dependency[site_id_0].insert(reaction_count);
                        site_reaction_dependency[site_id_1].insert(reaction_count);
                        reaction_count++;
                    }
                }
            }
        }
    }

    // // Resize initial propensities based on size of current reactions
    // initial_propensities.resize(initial_reactions.size());
    //
    // // initializing initial_propensities
    // for (unsigned int reaction_id = 0; reaction_id < initial_reactions.size(); reaction_id++) {
    //     initial_propensities[reaction_id] = compute_propensity(std::ref(initial_state), reaction_id);
    // }

}

void NanoParticle::compute_distance_matrix() {
    distance_matrix.resize(sites.size());
    for (unsigned int i = 0; i < sites.size(); i++) {
        distance_matrix[i].resize(sites.size());
        for (unsigned int j = 0; j < sites.size(); j++) {
            distance_matrix[i][j] = std::sqrt(site_distance_squared(sites[i], sites[j]));
        }
    }
}

std::vector<std::vector<int>> NanoParticle::compute_site_neighbors() {
    double threshold = interaction_radius_bound * interaction_radius_bound;
    std::vector<std::vector<int>> site_neighbors;
    site_neighbors.resize(sites.size());

    std::vector<int> buffer (sites.size());

    for (unsigned int i = 0; i < sites.size(); i++) {

        int count = 0;
        for (unsigned int j = 0; j < sites.size(); j++) {
            if ( i != j && // neighbors need to be distinct
                 site_distance_squared(sites[i], sites[j]) < threshold ) {
                buffer[j] = true;
                count += 1;
            }
            else {
                buffer[j] = false;
            }
        }

        site_neighbors[i].resize(count);

        count = 0;
        for (unsigned int j = 0; j < sites.size(); j++) {
            if (buffer[j]) {
                site_neighbors[i][count] = j;
                count += 1;
            }
        }
    }

    return site_neighbors;
}

void NanoParticle::update_state(
    std::vector<int> &state,
    Reaction reaction) {

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
}

void NanoParticle::update_reactions(
    std::vector<int> &state,
    std::vector<std::set<int>> &current_site_reaction_dependency,
    std::vector<Reaction> &current_reactions,
    Reaction reaction) {

    std::cerr << "---------\n"
              << "New "
              << reaction.interaction.number_of_sites
              << "-site reaction (involving sites "
              << reaction.site_id[0]
              << " and "
              << reaction.site_id[1]
              << ") fired: "
              << reaction.interaction.interaction_id
              << " - rate = "
              << reaction.rate
              << "\n";

    // Compute the new reactions based on the new states
    std::vector<Reaction> new_reactions;
    for ( int k = 0; k < reaction.interaction.number_of_sites; k++) {
        int site_id_0 = reaction.site_id[k];
        int site_0_state = reaction.interaction.right_state[k];
        int site_0_species_id = sites[site_id_0].species_id;

        // Add one site interactions
        std::vector<Interaction> available_interactions = one_site_interactions_map[site_0_species_id][site_0_state];
        for (unsigned int i = 0; i < available_interactions.size(); i++) {
            Interaction interaction = available_interactions[i];
            Reaction new_reaction = Reaction {
                                    .site_id = { (int) site_id_0, -1},
                                    .interaction = interaction,
                                    .rate = interaction.rate * one_site_interaction_factor
                                    };
            new_reactions.push_back(new_reaction);
        }

        // Add two site interactions
        for (unsigned int site_id_1 = 0; site_id_1 < sites.size(); site_id_1++) {
            if (site_id_0 != (int) site_id_1) {
                int site_1_state = state[site_id_1];
                int site_1_species_id = sites[site_id_1].species_id;
                std::vector<Interaction> available_interactions = two_site_interactions_map[site_0_species_id][site_1_species_id][site_0_state][site_1_state];
                double distance = distance_matrix[site_id_0][site_id_1];
                if (distance < interaction_radius_bound) {
                    for (unsigned int i = 0; i < available_interactions.size(); i++){
                        Interaction interaction = available_interactions[i];
                        Reaction new_reaction = Reaction {
                                                .site_id = { (int) site_id_0, (int) site_id_1 },
                                                .interaction = interaction,
                                                .rate = distance_factor_function(distance) * interaction.rate * two_site_interaction_factor
                                                };
                        new_reactions.push_back(new_reaction);
                    }
                }
            }
        }
    }
    // Find indexes of reactions to be removed
    std::set<int> reactions_to_remove;
    std::set<int>::iterator site_reaction_dependency_itr;
    std::set<int> reaction_dependency;
    for ( int k = 0; k < reaction.interaction.number_of_sites; k++) {
        std::cerr << "Adding to the remove list: ";
        reaction_dependency = current_site_reaction_dependency[reaction.site_id[k]];
        for (site_reaction_dependency_itr = reaction_dependency.begin();
             site_reaction_dependency_itr != reaction_dependency.end();
             site_reaction_dependency_itr++) {
              int reaction_id_to_remove = *site_reaction_dependency_itr;
              reactions_to_remove.insert(reaction_id_to_remove);

              std::cerr << *site_reaction_dependency_itr
                        << ", ";

              // Need to remove this interaction from the second site if this is a two_site interaction
              Reaction reaction_to_remove = current_reactions[reaction_id_to_remove];
              if (reaction_to_remove.interaction.number_of_sites == 2) {
                  current_site_reaction_dependency[reaction_to_remove.site_id[1]].erase(reaction_id_to_remove);
              }

        }

        std::cerr << "\n";
        current_site_reaction_dependency[reaction.site_id[k]].clear();
    }

    std::cerr << "There are currently "
              << current_reactions.size()
              << " reactions.\n"
              << reactions_to_remove.size()
              << " to be removed, and "
              << new_reactions.size()
              << " to be added.\n";

    std::cerr << "Removing: ";
    std::set<int>::iterator itr;
    for (itr = reactions_to_remove.begin(); itr != reactions_to_remove.end(); itr++) {
        std::cerr << *itr << ", ";
    }
    std::cerr << "\n";

    // Replace the reactions in current_reactions
    int reactions_replaced = 0;
    int n_reactions_to_remove = reactions_to_remove.size();
    itr = reactions_to_remove.begin();
    for (unsigned int i = 0; i < new_reactions.size(); i++) {
        Reaction new_reaction = new_reactions[i];
        if ( i >= n_reactions_to_remove){
            // If the number of new reactions is larger than the number of reactions to remove
            current_reactions.push_back(new_reaction);
            for (int k = 0; k < new_reaction.interaction.number_of_sites; k++) {
                current_site_reaction_dependency[new_reaction.site_id[k]].insert(current_reactions.size()-1);
            }
            std::cerr << "New interaction "
                      << current_reactions.size()-1
                      << ", involving sites "
                      << new_reaction.site_id[0]
                      << " and "
                      << new_reaction.site_id[1]
                      << "\n";
        } else {
            current_reactions[*itr] = new_reaction;
            for (int k = 0; k < new_reaction.interaction.number_of_sites; k++) {
                current_site_reaction_dependency[new_reaction.site_id[k]].insert(*itr);
            }
            std::cerr << "Replaced interaction "
                      << *itr
                      << "\n";
            reactions_replaced++;
            reactions_to_remove.erase(itr++);
            // std::advance(itr, 1);
        }
    }
    // reactions_to_remove.erase(reactions_to_remove.begin(), reactions_to_remove.begin()+reactions_replaced);

    std::cerr << "Removing: ";
    for (itr = reactions_to_remove.begin(); itr != reactions_to_remove.end(); itr++) {
        std::cerr << *itr << ", ";
    }
    std::cerr << "\n";

    // std::cerr << "Site Dependency map for site 426: ";
    // for (site_reaction_dependency_itr = current_site_reaction_dependency[426].begin();
    //      site_reaction_dependency_itr != current_site_reaction_dependency[426].end();
    //      site_reaction_dependency_itr++) {
    //     std::cerr << *site_reaction_dependency_itr
    //               << ", ";
    // }
    // std::cerr << "\n";

    std::cerr << "Num Reactions to remove: " << reactions_to_remove.size() << "\n";

    if (reactions_to_remove.size() > 0){
        int n_reactions_to_move = reactions_to_remove.size();
        int reactions_moved = 0;
        int reaction_id_to_move = (int) current_reactions.size()-1;

        std::set<int>::iterator reactions_moved_itr = reactions_to_remove.begin();
        while (reactions_moved < n_reactions_to_move) {
            // Could simplify this by enforcing order in current_site_reaction_dependency or by specifying a range
            auto result = reactions_to_remove.find(reaction_id_to_move);
            if (result != reactions_to_remove.end()){
                // Reaction to be moved is going to be removed, no point in removing it
                reaction_id_to_move--;
            } else if (reaction_id_to_move < *reactions_moved_itr) {
                break;
            } else {
                Reaction reaction_to_move = current_reactions[reaction_id_to_move];
                current_reactions[*reactions_moved_itr] = reaction_to_move;
                std::cerr << "Moving reaction "
                          << reaction_id_to_move
                          << " to "
                          << *reactions_moved_itr
                          << "\n";

                // find the reaction that was moved in the site reaction dependency vector and remap it
                for (int k = 0; k < reaction_to_move.interaction.number_of_sites; k++) {
                    int site_id = reaction_to_move.site_id[k];
                    // This number should always be in the list, else something has gone wrong
                    auto site_reaction_index_search = current_site_reaction_dependency[site_id].find(reaction_id_to_move);

                    if (site_reaction_index_search != current_site_reaction_dependency[site_id].end()){
                        current_site_reaction_dependency[site_id].erase(reaction_id_to_move);
                        // int site_reaction_index = site_reaction_index_search - current_site_reaction_dependency[site_id].begin();

                        std::cerr << "Remapping reaction "
                                  << *site_reaction_index_search
                                  << " to "
                                  << *reactions_moved_itr
                                  << "\n";
                        current_site_reaction_dependency[site_id].insert(*reactions_moved_itr);
                    } else {
                      // Throw some sort of error
                      std::cerr << "Could not find reaction "
                                << reaction_id_to_move
                                << " in the site reaction dependency map for site "
                                << site_id
                                << "\n";

                      for (itr = current_site_reaction_dependency[site_id].begin();
                           itr != current_site_reaction_dependency[site_id].end();
                           itr++) {
                             std::cerr << *itr
                                       << ", ";
                      }
                      std::cerr << "\n";
                      raise(SIGINT);
                    }
                }
                reaction_id_to_move--;
                reactions_moved_itr++;
                reactions_moved++;
            }
        }
        // for (unsigned int i = new_reactions.size(); i < reactions_to_remove.size(); i++){
        //     // take a reaction from the end of the list and put it into the place of the reaction to remove
        //     int reaction_index_to_move = current_reactions.size()-i;
        //     Reaction reaction_to_move = current_reactions[reaction_index_to_move];
        //     current_reactions[reactions_to_remove[i]] = reaction_to_move;
        //
        //     // find the reaction that was moved in the site reaction dependency vector and remap it
        //     for (int k = 0; k < reaction_to_move.interaction.number_of_sites; k++) {
        //         int site_id = reaction_to_move.site_id[k];
        //         // This number should always be in the list, else something has gone wrong
        //         auto site_reaction_index_search = std::find(current_site_reaction_dependency[site_id].begin(), current_site_reaction_dependency[site_id].end(), reaction_index_to_move);
        //
        //         if (site_reaction_index_search != current_site_reaction_dependency[site_id].end()){
        //             int site_reaction_index = site_reaction_index_search - current_site_reaction_dependency[site_id].begin();
        //             current_site_reaction_dependency[site_id][site_reaction_index] = reactions_to_remove[i];
        //         } else {
        //           // Throw some sort of error
        //           std::cerr << "Could not find reaction "
        //                     << reaction_index_to_move
        //                     << " in the site reaction dependency map for site "
        //                     << site_id
        //                     << "\n";
        //
        //           raise(SIGINT);
        //         }
        //     }
        // }
        current_reactions.resize(current_reactions.size() - n_reactions_to_move);
    }
}

TrajectoriesSql NanoParticle::history_element_to_sql(
    int seed,
    HistoryElement history_element) {

    Reaction reaction = history_element.reaction;
    return TrajectoriesSql {
        .seed = seed,
        .step = history_element.step,
        .time = history_element.time,
        .site_id_1 = reaction.site_id[0],
        .site_id_2 = reaction.site_id[1],
        .interaction_id = reaction.interaction.interaction_id
    };
}

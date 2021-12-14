#pragma once
#include "../core/sql.h"
#include "../core/simulation.h"
#include "sql_types.h"
#include <vector>
#include <cmath>
#include <functional>

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

struct Interaction {
    // either 1 site or two site interaction
    int number_of_sites;
    int species_id[2];
    int left_state[2];
    int right_state[2];
    double rate;
};

// For this nano particle simulator, a reaction consists
// of upto two sites and the interaction id.
// if it is an internal interaction, site_id_2 will be -1
// In a reaction, the sites must be within the interaction radius bound.

struct Reaction {
    int site_id[2];
    int interaction_id;
    double distance;
};

struct NanoParticleParameters {};

struct NanoParticle {
    // maps a species index to the number of degrees of freedom
    std::vector<int> degrees_of_freedom;

    // maps site index to site data
    std::vector<Site> sites;

    // maps site ids to reaction ids involving the site.
    std::vector<std::vector<int>> site_reaction_dependency;

    // maps interaction index to interaction data
    std::vector<Interaction> interactions;

    // initial state of the simulations.
    // initial_state[i] is a local degree of freedom
    // from the species at site i.
    std::vector<int> initial_state;

    std::vector<double> initial_propensities;

    // list mapping reaction_ids to reactions
    std::vector<Reaction> reactions;

    double one_site_interaction_factor;
    double two_site_interaction_factor;
    double interaction_radius_bound;

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

    double compute_propensity(
        std::vector<int> &state,
        int reaction_id);

    void update_state(
        std::vector<int> &state,
        int reaction_id);

    void update_propensities(
        std::function<void(Update update)> update_function,
        std::vector<int> &state,
        int next_reaction_id
        );

    // convert a history element as found a simulation to history
    // to a SQL type.
    TrajectoriesSql history_element_to_sql(
        int seed,
        int step,
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
    interactions.resize(metadata_row.number_of_interactions);
    while(std::optional<InteractionSql> maybe_interaction_row =
          interactions_reader.next()) {

        InteractionSql interaction_row = maybe_interaction_row.value();
        interactions[interaction_row.interaction_id] = {
            .number_of_sites = interaction_row.number_of_sites,
            .species_id      = { interaction_row.species_id_1, interaction_row.species_id_2},
            .left_state      = { interaction_row.left_state_1, interaction_row.left_state_2},
            .right_state     = { interaction_row.right_state_1, interaction_row.right_state_2},
            .rate            = interaction_row.rate
        };
    }

    // initialize initial_state
    initial_state.resize(metadata_row.number_of_sites);

    while(std::optional<InitialStateSql> maybe_initial_state_row =
          initial_state_reader.next()) {
        InitialStateSql initial_state_row = maybe_initial_state_row.value();
        initial_state[initial_state_row.site_id] = initial_state_row.degree_of_freedom;
    }

    compute_reactions();
    initial_propensities.resize(reactions.size());

    // initializing initial_propensities
    for (unsigned int reaction_id = 0; reaction_id < reactions.size(); reaction_id++) {
        initial_propensities[reaction_id] = compute_propensity(std::ref(initial_state), reaction_id);
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
            if (site_distance_squared(sites[i], sites[j]) < threshold) {
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

void NanoParticle::compute_reactions() {
    // compute all possible reactions and a mapping of site ids to the
    // reaction ids involving a site. Since we expect the number of distinct interactions
    // to be quite small, the number of reactions involving two fixed sites is quite small.
    // this means that we don't loose much by indexing dependency by site id rather than
    // reaction id, as the intersection will be small.
    int reaction_count = 0;
    std::vector<std::vector<int>> site_neighbors = compute_site_neighbors();
    std::vector<int> site_reaction_dependency_counter;
    site_reaction_dependency_counter.resize(sites.size());

    // counting one site interactions
    for ( unsigned int site_id = 0;
          site_id < sites.size();
          site_id++ ) {
        for ( unsigned int interaction_id = 0;
              interaction_id < interactions.size();
              interaction_id++ ) {
            int species_id = sites[site_id].species_id;
            if ((interactions[interaction_id].number_of_sites == 1) &&
                (species_id == interactions[interaction_id].species_id[0])) {
                reaction_count += 1;
                site_reaction_dependency_counter[site_id] += 1;
            }
        }
    }

    // counting two site interactions
    for ( unsigned int site_id_0 = 0;
          site_id_0 < sites.size();
          site_id_0++) {
        for ( unsigned int j = 0;
              j < site_neighbors[site_id_0].size();
              j++) {

            unsigned int site_id_1 = site_neighbors[site_id_0][j];

            for ( unsigned int interaction_id = 0;
                  interaction_id < interactions.size();
                  interaction_id++) {

                int species_id_0 = sites[site_id_0].species_id;
                int species_id_1 = sites[site_id_1].species_id;

                if ((interactions[interaction_id].number_of_sites == 2) &&
                    (interactions[interaction_id].species_id[0] == species_id_0) &&
                    (interactions[interaction_id].species_id[1] == species_id_1)) {

                    reaction_count += 1;
                    site_reaction_dependency_counter[site_id_0] += 1;
                    site_reaction_dependency_counter[site_id_1] += 1;
                }
            }
        }
    }


    reactions.resize(reaction_count);

    // reseting reaction_count and site_reaction_dependency_count
    // and resizing the entries in site_reaction_dependency
    reaction_count = 0;
    for (unsigned int i = 0; i < site_reaction_dependency.size(); i++) {
        site_reaction_dependency[i].resize(site_reaction_dependency_counter[i]);
        site_reaction_dependency_counter[i] = 0;
    }

    // setting one site reactions
    for ( unsigned int site_id = 0;
          site_id < sites.size();
          site_id++) {
        for ( unsigned int interaction_id = 0;
              interaction_id < interactions.size();
              interaction_id++) {
            int species_id = sites[site_id].species_id;
            if ((interactions[interaction_id].number_of_sites == 1) &&
                (species_id == interactions[interaction_id].species_id[0])) {

                reactions[reaction_count] = {
                    .site_id = { (int) site_id, -1},
                    .interaction_id = (int) interaction_id,
                    .distance = 0
                };
                site_reaction_dependency[site_id][
                    site_reaction_dependency_counter[site_id]] = reaction_count;

                reaction_count += 1;
                site_reaction_dependency_counter[site_id] += 1;

            }
        }
    }

    // setting two site interactions
    for ( unsigned int site_id_0 = 0;
          site_id_0 < sites.size();
          site_id_0++) {
        for ( unsigned int j = 0;
              j < site_neighbors[site_id_0].size();
              j++) {

            unsigned int site_id_1 = site_neighbors[site_id_0][j];

            for ( unsigned int interaction_id = 0;
                  interaction_id < interactions.size();
                  interaction_id++) {

                Site site_0 = sites[site_id_0];
                Site site_1 = sites[site_id_1];
                int species_id_0 = site_0.species_id;
                int species_id_1 = site_1.species_id;


                if ((interactions[interaction_id].number_of_sites == 2) &&
                    (interactions[interaction_id].species_id[0] == species_id_0) &&
                    (interactions[interaction_id].species_id[1] == species_id_1)) {

                    reactions[reaction_count] = {
                        .site_id = { (int) site_id_0, (int) site_id_1 },
                        .interaction_id = (int) interaction_id,
                        .distance = std::sqrt(site_distance_squared(site_0, site_1))

                    };

                    site_reaction_dependency[site_id_0][
                        site_reaction_dependency_counter[site_id_0]] = reaction_count;

                    site_reaction_dependency[site_id_1][
                        site_reaction_dependency_counter[site_id_1]] = reaction_count;

                    reaction_count += 1;
                    site_reaction_dependency_counter[site_id_0] += 1;
                    site_reaction_dependency_counter[site_id_1] += 1;
                }
            }
        }
    }
}

double NanoParticle::compute_propensity(
    std::vector<int> &state,
    int reaction_id) {

    Interaction interaction = interactions[reactions[reaction_id].interaction_id];

    if ( interaction.number_of_sites == 1) {
        // internal interaction
        int site_id = reactions[reaction_id].site_id[0];
        if (interaction.left_state[0] == state[site_id]) {
            return (
                interaction.rate *
                one_site_interaction_factor);
        }
        else {
            return 0;
        }


    }
    else {
        // two site interaction
        int site_id_0 = reactions[reaction_id].site_id[0];
        int site_id_1 = reactions[reaction_id].site_id[1];

        if ( (interaction.left_state[0] == state[site_id_0]) &&
             (interaction.left_state[1] == state[site_id_1]) ) {

            // TODO: add distance as an attribute to the reaction struct
            double distance = reactions[reaction_id].distance;

            // by the way we constructed the reactions, diff will always be >= 0
            double distance_factor = 1 - ( distance / interaction_radius_bound );
            return (
                interaction.rate *
                distance_factor *
                two_site_interaction_factor);


        }
        else {
            return 0;
        }
    }
}


void NanoParticle::update_state(
    std::vector<int> &state,
    int reaction_id) {
    Reaction reaction = reactions[reaction_id];

    Interaction interaction = interactions[
        reaction.interaction_id];

    for (int k = 0; k < interaction.number_of_sites; k++) {
        state[reaction.site_id[k]] = interaction.right_state[k];
    }
}


void NanoParticle::update_propensities(
    std::function<void(Update update)> update_function,
    std::vector<int> &state,
    int next_reaction_id
    ) {
    Reaction reaction = reactions[next_reaction_id];
    Interaction interaction = interactions[
        reaction.interaction_id];

    for ( int k = 0; k < interaction.number_of_sites; k++) {
        for ( unsigned int i = 0;
              i < site_reaction_dependency[reaction.site_id[k]].size();
              i++ ) {
            int reaction_id = site_reaction_dependency[reaction.site_id[k]][i];
            double new_propensity = compute_propensity(std::ref(state), reaction_id);


            update_function( Update {
                    .index = (unsigned long int) reaction_id,
                    .propensity = new_propensity});
        }
    }
}


TrajectoriesSql NanoParticle::history_element_to_sql(
    int seed,
    int step,
    HistoryElement history_element) {

    Reaction reaction = reactions[history_element.reaction_id];
    return TrajectoriesSql {
        .seed = seed,
        .step = step,
        .time = history_element.time,
        .site_id_1 = reaction.site_id[0],
        .site_id_2 = reaction.site_id[1],
        .interaction_id = reaction.interaction_id
    };
}


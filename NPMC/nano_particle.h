#pragma once
#include "../core/sql.h"
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
    int species_id_1;
    int species_id_2;
    int left_state_1;
    int left_state_2;
    int right_state_1;
    int right_state_2;
    double rate;
};

// For this nano particle simulator, a reaction consists
// of upto two sites and the interaction id.
// if it is an internal interaction, site_id_2 will be -1
// In a reaction, the sites must be within the interaction radius bound.

struct Reaction {
    int site_id_1;
    int site_id_2;
    int interaction_id;
};

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
        SqlConnection &initial_state_database);

    // maps a site index to the indices of its neighbors
    // within the spatial decay radius
    std::vector<std::vector<int>> compute_site_neighbors();

    void compute_reactions();

    double compute_propensity(
        std::vector<int> &state,
        int reaction_id);

};


NanoParticle::NanoParticle(
    SqlConnection &nano_particle_database,
    SqlConnection &initial_state_database) {

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
        // computing the type of interaction
        int number_of_sites = 0;
        if (interaction_row.left_state_2 < 0)
            number_of_sites = 1;
        else
            number_of_sites = 2;

        interactions[interaction_row.interaction_id] = {
            .number_of_sites = number_of_sites,
            .species_id_1    = interaction_row.species_id_1,
            .species_id_2    = interaction_row.species_id_2,
            .left_state_1    = interaction_row.left_state_1,
            .left_state_2    = interaction_row.left_state_2,
            .right_state_1   = interaction_row.right_state_1,
            .right_state_2   = interaction_row.right_state_2,
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
                (species_id == interactions[interaction_id].species_id_1)) {
                reaction_count += 1;
                site_reaction_dependency_counter[site_id] += 1;
            }
        }
    }

    // counting two site interactions
    for ( unsigned int site_id_1 = 0;
          site_id_1 < sites.size();
          site_id_1++) {
        for ( unsigned int j = 0;
              j < site_neighbors[site_id_1].size();
              j++) {

            unsigned int site_id_2 = site_neighbors[site_id_1][j];

            for ( unsigned int interaction_id = 0;
                  interaction_id < interactions.size();
                  interaction_id++) {

                int species_id_1 = sites[site_id_1].species_id;
                int species_id_2 = sites[site_id_2].species_id;

                if ((interactions[interaction_id].number_of_sites == 2) &&
                    (interactions[interaction_id].species_id_1 == species_id_1) &&
                    (interactions[interaction_id].species_id_2 == species_id_2)) {

                    reaction_count += 1;
                    site_reaction_dependency_counter[site_id_1] += 1;
                    site_reaction_dependency_counter[site_id_2] += 1;
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
                (species_id == interactions[interaction_id].species_id_1)) {

                reactions[reaction_count] = {
                    .site_id_1 = (int) site_id,
                    .site_id_2 = -1,
                    .interaction_id = (int) interaction_id
                };
                site_reaction_dependency[site_id][
                    site_reaction_dependency_counter[site_id]] = reaction_count;

                reaction_count += 1;
                site_reaction_dependency_counter[site_id] += 1;

            }
        }
    }

    // setting two site interactions
    for ( unsigned int site_id_1 = 0;
          site_id_1 < sites.size();
          site_id_1++) {
        for ( unsigned int j = 0;
              j < site_neighbors[site_id_1].size();
              j++) {

            unsigned int site_id_2 = site_neighbors[site_id_1][j];

            for ( unsigned int interaction_id = 0;
                  interaction_id < interactions.size();
                  interaction_id++) {

                int species_id_1 = sites[site_id_1].species_id;
                int species_id_2 = sites[site_id_2].species_id;

                if ((interactions[interaction_id].number_of_sites == 2) &&
                    (interactions[interaction_id].species_id_1 == species_id_1) &&
                    (interactions[interaction_id].species_id_2 == species_id_2)) {

                    reactions[reaction_count] = {
                        .site_id_1 = (int) site_id_1,
                        .site_id_2 = (int) site_id_2,
                        .interaction_id = (int) interaction_id
                    };

                    site_reaction_dependency[site_id_1][
                        site_reaction_dependency_counter[site_id_1]] = reaction_count;

                    site_reaction_dependency[site_id_2][
                        site_reaction_dependency_counter[site_id_2]] = reaction_count;

                    reaction_count += 1;
                    site_reaction_dependency_counter[site_id_1] += 1;
                    site_reaction_dependency_counter[site_id_2] += 1;

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
        int site_id = reactions[reaction_id].site_id_1;
        if (interaction.left_state_1 == state[site_id]) {
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
        int site_id_1 = reactions[reaction_id].site_id_1;
        int site_id_2 = reactions[reaction_id].site_id_2;

        if ( (interaction.left_state_1 == state[site_id_1]) &&
             (interaction.left_state_2 == state[site_id_2]) ) {
            Site site_1 = sites[site_id_1];
            Site site_2 = sites[site_id_2];
            double distance = std::sqrt(site_distance_squared(site_1, site_2));


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

#pragma once
#include "../core/sql.h"
#include "sql_types.h"
#include <vector>



struct Site {
    double x;
    double y;
    double z;
    unsigned int species_id;
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

struct NanoParticle {
    // maps a species index to the number of degrees of freedom
    std::vector<int> degrees_of_freedom;

    // maps site index to site data
    std::vector<Site> sites;

    // maps a site index to the indices of its neighbors
    // within the spatial decay radius
    std::vector<std::vector<int>> site_neighbors;

    // maps interaction index to interaction data
    std::vector<Interaction> interactions;

    double one_site_interaction_factor;
    double two_site_interaction_factor;
    double spatial_decay_radius;

    // constructor
    NanoParticle(
        SqlConnection &nano_particle_database,
        SqlConnection &initial_state_database);

    void compute_site_neighbors();
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

    // sql readers
    SqlReader<SpeciesSql> species_reader(species_statement);
    SqlReader<SiteSql> site_reader(site_statement);
    SqlReader<InteractionSql> interactions_reader(interactions_statement);
    SqlReader<MetadataSql> metadata_reader(metadata_statement);
    SqlReader<FactorsSql> factors_reader(factors_statement);

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
    spatial_decay_radius = factor_row.spatial_decay_radius;


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
    site_neighbors.resize(metadata_row.number_of_sites);

    while(std::optional<SiteSql> maybe_site_row =
          site_reader.next()) {

        SiteSql site_row = maybe_site_row.value();
        sites[site_row.site_id] = {
            .x = site_row.x,
            .y = site_row.y,
            .z = site_row.z,
            .species_id = site_row.species_id };
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

    compute_site_neighbors();

}

void NanoParticle::compute_site_neighbors() {
    // For all sites, compute the sites which are within spatial decay radius
    // For this kind of computation, there is always a trade off.
    // do you allocate all the arrays in one chunk, which means allocating more mem
    // overall, or do you do a bunch of individual allocations. Since we run this once
    // at the start and want the mem footprint to be as small as possible, we
    // go for the second option.
    double threshold = spatial_decay_radius * spatial_decay_radius;

    std::vector<unsigned int> buffer (sites.size());

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
}

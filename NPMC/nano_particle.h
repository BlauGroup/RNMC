#pragma once
#include "../core/sql.h"
#include <vector>


struct Site {
    double x;
    double y;
    double z;
    int species_id;
};

struct NanoParticle {
    std::vector<Site> sites;

    NanoParticle(
        SqlConnection *nano_particle_database,
        SqlConnection *initial_state_database);
};


NanoParticle::NanoParticle(
    SqlConnection *nano_particle_database,
    SqlConnection *initial_state_database) {
};

/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#ifndef RNMC_NPMC_TYPES_H
#define RNMC_NPMC_TYPES_H

struct NanoParticleParameters
{
    bool isCheckpoint;
};

struct Interaction
{
    // either 1 site or two site interaction
    int interaction_id;
    int number_of_sites;
    int species_id[2];
    int left_state[2];
    int right_state[2];

    // the units of rate depend on number_of_sites. If number_of_sites = 1, then
    // rate has units 1 / s. If number of sites = 2, then rate has units 1 / s m^6.
    double rate;
};

struct NanoSite
{
    double x;
    double y;
    double z;
    int species_id;
};

// For this nano particle simulator, a reaction consists
// of upto two sites and the interaction id.
// if it is an internal interaction, site_id_2 will be -1
// In a reaction, the sites must be within the interaction radius bound.
struct NanoReaction
{
    int site_id[2];
    Interaction interaction;

    // rate has units 1 / s
    double rate;
};

#endif
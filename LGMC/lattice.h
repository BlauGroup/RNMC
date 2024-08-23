/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.

Lattice from https://spparks.sandia.gov/
---------------------------------------------------------------------- */

#ifndef RNMC_LATTICE_H
#define RNMC_LATTICE_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <map>
#include <cassert>
#include <iostream>

#include "stdio.h"
#include "stdlib.h"
#include <string>
#include "math.h"

#define SPECIES_EMPTY 0
#define DELTALOCAL 10000
#define DELTA 32768
#define EPSILON 0.0001

#define FLERR __FILE__, __LINE__

#ifndef MIN
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

struct Site
{

    Site() {}

    Site(uint32_t i_in, uint32_t j_in, // custom constructor for convenience
         uint32_t k_in,
         int species_in,
         bool can_adsorb_in) : i(i_in), j(j_in), k(k_in),
                               species(species_in), can_adsorb(can_adsorb_in)
    {
    }

    uint32_t i;      // site location on lattice
    uint32_t j;
    uint32_t k;
    int species;     // species at site
    bool can_adsorb; // is the site in contact with the electrolyte?
};

class Lattice
{
private:
    /* ----------------------- handle memory ---------------------------------- */

    void *smalloc(int nbytes, const char *name); // safe allocate

    template <typename TYPE>
    TYPE *create(TYPE *&array, int n, const char *name); // create 1D array

    template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2, const char *name); // create 2D array

    template <typename TYPE>
    void destroy(TYPE **array); // destroy 2D array

    void sfree(void *ptr); // safe free

    /* ----------------------- structural information ------------------------ */

    int ilo, ihi, klo, khi, jlo, jhi; // geometry info neighbors
                                      // 0 = non-periodic, 1 = periodic

    int nsites;                       // number of sites
    int nmax;                         // max # sites, idneigh, numneigh

    int maxneigh;                     // max neigh per site
    float maxz;                       // the largest dist of lattice 

public:
    float latconst;                   // lattice constant
    bool is_xperiodic, is_yperiodic, is_zperiodic;
    float xlo, xhi, ylo,              // bounding of box
        yhi, zlo, zhi;                // (yscale * read in value)

    std::unordered_map<int, Site> sites;         // list of Sites for lattice
    std::unordered_map<int, uint32_t *> idneigh; // neighbor IDs for each site
    std::unordered_map<int, uint32_t> numneigh;  // # of neighbors of each site
    std::unordered_map<int, char> edges;
    std::map<std::tuple<uint32_t, uint32_t, uint32_t>, int> loc_map; // Mapping from site location (i,j,k) to site ID
    bool isCheckpoint;

    Lattice(float latconst_in);

    Lattice(float latconst_in, int ihi_in,
            int jhi_in, int khi_in);

    Lattice(const Lattice &other); // copy constructor

    ~Lattice();

    void structured_lattice();

    void structured_connectivity();

    void offsets_3d(int **cmapone);

    void add_site(uint32_t i_in, uint32_t j_in, uint32_t k_in,
                  bool can_adsorb_in, bool update_neighbors_in, bool meta_neighbors_in);

    void delete_site(int id);

    void update_neighbors(uint32_t n, bool meta_neighbors_in);

    float get_latconst();

    float get_maxz();
};

#endif
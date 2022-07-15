#ifndef LATTICE_H
#define LATTICE_H

#include "memory.h"
#include <vector>
#include <unordered_map>
#include <string>

namespace LGMC_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif
#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

struct Site {

    Site() {}                               // default constructor does nothing

    Site(uint32_t i_in, uint32_t j_in,      // custom constructor for convenience 
        uint32_t k_in, float x_in, float y_in, float z_in,
        int species_in,
        bool can_adsorb_in) :
        i(i_in), j(j_in), k(k_in),
        x(x_in), y(y_in), z(z_in), 
        species(species_in), can_adsorb(can_adsorb_in) { }

    uint32_t i;                             // site location on lattice
    uint32_t j;
    uint32_t k;
    float x;                                // location in space
    float y;
    float z;
    int species;                            // species at site
    bool can_adsorb;                        // is the site in contact with the electrolyte?

};

class Lattice {

private:
    Memory *memory;                             // handle lattice memory
    
    /* ----------------------- structural information ------------------------ */
    
    float latconst;                               // lattice constant
    int boxxlo,boxxhi,boxylo,                   // bounding of box
    boxyhi,boxzlo,boxzhi;                       // (yscale * read in value)
    int xlo,xhi,ylo,yhi,zlo,zhi;                // geometry info neighbors
    bool is_xperiodic, is_yperiodic, is_zperiodic;          // 0 =   non-periodic, 1 = periodic
    
    int nsites;                                 // number of sites
    int nmax;                                   // max # of sites per-site arrays can store

    int maxneigh;                               // max neighbors per site

public: 

    std::vector<Site> sites;                                // list of Sites for lattice
    std::vector<uint32_t*> idneigh;                         // neighbor IDs for each site
    std::vector<uint32_t> numneigh;                         // # of neighbors of each site

    Lattice(float latconst_in, 
        int boxxlo_in, int boxxhi_in, int boxylo_in,
        int boxyhi_in, int boxzlo_in, int boxzhi_in, 
        bool xperiodic_in, bool yperiodic_in, bool zperiodic_in);   

    ~Lattice();

    void structured_lattice();
    
    void structured_connectivity();
    
    void offsets_3d(int **cmapone);
    
    void add_site(uint32_t i_in, uint32_t j_in, 
                  uint32_t k_in, float x_in, float y_in, float z_in,
                  bool can_adsorb_in);

    void update_neighbors(unit32_t n);
};

}

#endif

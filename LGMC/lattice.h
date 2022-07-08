#ifndef LATTICE_H
#define LATTICE_H

#include "memory.h"
#include <vector>
#include <unordered_map>
#include <string>
#include "../core/dispatcher.h"
#include "../GMC/reaction_network.h"
#include "../GMC/sql_types.h"

namespace LGMC_NS {

// universal defines inside namespace

#define FLERR __FILE__,__LINE__

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif
#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

class Lattice {

private:
    Memory *memory;                             // handle lattice memory
    
    /* ----------------------- structural information ------------------------ */
    
    int latconst;                               // lattice constant
    int boxxlo,boxxhi,boxylo,                   // bounding of box
    boxyhi,boxzlo,boxzhi;                       // (yscale * read in value)
    int xlo,xhi,ylo,yhi,zlo,zhi;                // geometry info neighbors
    uint32_t xperiodic,yperiodic,zperiodic;          // 0 = non-periodic, 1 = periodic
    
    int nsites;                                 // number of sites
    int nmax;                                   // max # of sites per-site arrays can store

    struct Site {

        Site() {}                               // default constructor does nothing

        
        Site(uint32_t x_in, uint32_t y_in,      // custom constructor for convenience 
            uint32_t z_in, int species_in) :
            x(x_in), y(y_in), z(z_in), 
            species(species_in) { }

        uint32_t x;                             // site location on lattice
        uint32_t y;
        uint32_t z;
        int species;                            // species at site

    };

    Site *sites;                                // list of Sites for lattice

    int maxneigh;                               // max neighbors per site
    uint32_t **idneigh;                         // neighbor IDs for each site
    uint32_t *numneigh;                         // # of neighbors of each site

    /* -------------------------- kMC information ------------------------------ */
    
    std::unordered_map<std::string,     // lattice propensities as site neighbor pair
                        std::vector<double> > props;          // -1 single site, -2 gillespie
                                                 

    std::vector<std::vector<int>>               // dependency list without products  
    lat_dependents;         

    int prop_sum;                               // running total of propensities
    
    <Solver, Model, Parameters, TrajectoriesSql>Dispatcher *dis_ptr;                        // pointer to gillespie dispatcher object
    
public: 
    Lattice(int latconst_in, <Solver, Model, Parameters, TrajectoriesSql>Dispatcher *ptr_in, 
        int boxxlo_in, int boxxhi_in, int boxylo_in,
        int boxyhi_in, int boxzlo_in, int boxzhi_in);   

    ~Lattice();

    void structured_lattice();
    
    void structured_connectivity();
    
    void offsets_3d(int **cmapone);
    
    void add_site(uint32_t n, uint32_t x, uint32_t y, uint32_t z);

    void grow(uint32_t n);
    
    void update_propensity(int site_one, int site_two, Lat_Reaction &reaction);

    void update(int site_one, int site_two);

    void relevant_react(int site);

    void clear_site(int site);
    
};

}

#endif

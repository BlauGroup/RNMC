#ifndef LATTICE_H
#define LATTICE_H

#include "memory.h"

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
    
    double latconst;                            // lattice constant
    double boxxlo,boxxhi,boxylo,                // bounding of box
    boxyhi,boxzlo,boxzhi;                       // (yscale * read in value)
    int xlo,xhi,ylo,yhi,zlo,zhi;                // geometry info neighbors
    int xperiodic,yperiodic,zperiodic;          // 0 = non-periodic, 1 = periodic
    
    int nsites;                                 // number of sites
    int nmax;                                   // max # of sites per-site arrays can store
    
    tagint *id;                                 // global ID of site
    double **xyz;                               // coords of site
    int *species;                               // ID of species at each lattice site

    
    int maxneigh;                               // max neighbors per site
    tagint **idneigh;                           // neighbor IDs for each site
    int *numneigh;                              // # of neighbors of each site

    /* -------------------------- kMC information ------------------------------ */
    int *state;                                 // count of all species
    int **propensities;                         // propensities of all site combinations
    
    
    
    
    
public: 
    Lattice(int, char **arg);               // style, dimension, latconst
    ~Lattice();

    void structured_lattice();
    
    void structured_connectivity();
    
    void offsets_3d(int **cmapone);
    
    void add_site(tagint n, double x, double y, double z);

    void grow(int n);
    
    //TODO: implement
    void init();
    
    void update();
    
    void compute_propensity();
    
    

};

}

#endif

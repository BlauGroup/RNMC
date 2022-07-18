#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <unordered_map>
#include <string>
#include <map>

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

    /* ----------------------- handle memory ---------------------------------- */
    
    void *smalloc(int nbytes, const char *name);              // safe allocate         

    template <typename TYPE>
    TYPE *create(TYPE *&array, int n, const char *name);         // create 1D array

    template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2, const char *name); // create 2D array

    template <typename TYPE>
    void destroy(TYPE **array);                                  // destroy 2D array

    void sfree(void *ptr);                                        // safe free
    
    /* ----------------------- structural information ------------------------ */
    
    float latconst;                               // lattice constant
    float boxxlo,boxxhi,boxylo,                   // bounding of box
    boxyhi,boxzlo,boxzhi;                       // (yscale * read in value)
    int xlo,xhi,ylo,yhi,zlo,zhi;                // geometry info neighbors
    bool is_xperiodic, is_yperiodic, is_zperiodic;          
    
    int nsites;                                 // number of sites
    int nmax;                                   // max # sites, idneigh, numneigh can store at a given time

    int maxneigh;                               // max neighbors per site

public: 

    std::vector<Site> sites;                                // list of Sites for lattice
    std::vector<uint32_t*> idneigh;                         // neighbor IDs for each site
    std::vector<uint32_t> numneigh;                         // # of neighbors of each site
    std::vector<uint32_t> edge;                             // ids of sites on the edge that can adsorb

    std::map<std::tuple<uint32_t, uint32_t, uint32_t>, int> loc_map;  // Mapping from site location to site ID

    Lattice(float latconst_in, 
        float boxxlo_in, float boxxhi_in, float boxylo_in,
        float boxyhi_in, float boxzlo_in, float boxzhi_in, 
        bool xperiodic_in, bool yperiodic_in, bool zperiodic_in);   
    
    Lattice(const Lattice& other);                          // copy constructor

    ~Lattice();

    void structured_lattice();
    
    void structured_connectivity();
    
    void offsets_3d(int **cmapone);
    
    void add_site(uint32_t i_in, uint32_t j_in, 
                  uint32_t k_in, float x_in, float y_in, float z_in,
                  bool can_adsorb_in, bool update_neighbors_in, bool meta_neighbors_in);

    void update_neighbors(uint32_t n, bool meta_neighbors_in);

    
};

}

#endif

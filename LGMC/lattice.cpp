#include "lattice.h"

Lattice::Lattice(float latconst_in) {

    isCheckpoint = true;
    latconst = latconst_in;
    
    // region of simulation input * lattice spacing
    xlo = 0;
    xhi = 1 * latconst;
    ylo = 0;
    yhi = 1 * latconst;
    zlo = 0;
    zhi = 1 * latconst;
    
    // 0 = non-periodic, 1 = periodic
    is_xperiodic = true;
    is_yperiodic = true;
    is_zperiodic = false;
    
    nsites = 0;
    maxz = 0;
    nmax = DELTA;

    maxneigh = 6;
        
} // Lattice()

/* ---------------------------------------------------------------------- */

Lattice::Lattice(float latconst_in, int ihi_in, 
        int jhi_in, int khi_in)  {
    
    isCheckpoint = false;
    latconst = latconst_in;
    
    // region of simulation input * lattice spacing
    xlo = 0;
    xhi = ihi_in * latconst;
    ylo = 0;
    yhi = jhi_in * latconst;
    zlo = 0;
    zhi = khi_in * latconst;
    
    // 0 = non-periodic, 1 = periodic
    is_xperiodic = true;
    is_yperiodic = true;
    is_zperiodic = false;
    
    nsites = 0;
    maxz = 0;
    nmax = DELTA;

    maxneigh = 6;

    // create sites on lattice
    structured_lattice();
    
    // set neighbors of each site
    structured_connectivity();

} // Lattice()

/* ---------------------------------------------------------------------- */

Lattice::Lattice(const Lattice& other) {

    latconst = other.latconst;                               
    ilo = other.ilo;
    ihi = other.ihi;
    jlo = other.jlo;               
    jhi = other.jhi;
    klo = other.klo;
    khi = other.khi;                       
    xlo = other.xlo;
    xhi = other.xhi;
    ylo = other.ylo;
    yhi = other.yhi;
    zlo = other.zlo;
    zhi = other.zhi;                
    is_xperiodic = other.is_xperiodic;
    is_yperiodic = other.is_yperiodic;
    is_zperiodic = other.is_zperiodic;         
    
    for(auto it = other.sites.begin(); it != other.sites.end(); it++) {
        sites[it->first] = it->second;
    }

    nsites = other.nsites;                               
    nmax = other.nmax;                                

    maxneigh = other.maxneigh;                                             
        
    for(auto it = other.numneigh.begin(); it != other.numneigh.end(); it++) {
        numneigh[it->first] = it->second;
    }

    for(auto it = other.edges.begin(); it != other.edges.end(); it++) {
        edges[it->first] = it->second;
    }                       

    for(int i = 0; i < int(numneigh.size()); i++) {
        uint32_t* neighi;
        create(neighi, maxneigh, "create:neighi");
 
        for(size_t j = 0; j < other.numneigh.at(i); j++) {
            neighi[j] = other.idneigh.at(i)[j];
        }
        idneigh[i] = neighi;
    }                           


} // Lattice()

/* ---------------------------------------------------------------------- */

Lattice::~Lattice() {

    for(size_t i = 0; i < idneigh.size(); i++) {
        sfree(idneigh[i]);
    }
} // ~Lattice()


/* ---------------------------------------------------------------------- */

void Lattice::structured_lattice() {
    
    // if not fully periodic IDs may be non-contiguous and/or ordered irregularly
    uint32_t nx, ny, nz;
    nx = (xhi - xlo / latconst);
    ny = (yhi - ylo / latconst);
    nz = (zhi - zlo / latconst);

    // if dim is periodic:
    //    lattice origin = lower box boundary
    //    loop bounds = 0 to N-1
    // if dim is non-periodic:
    //   lattice origin = 0.0
    //   loop bounds = enough to tile box completely, with all basis atoms
    
    if (is_xperiodic) {
        ilo = 0;
        ihi = nx-1;
      }
    else {
        ilo = (xlo / latconst);
        while ((ilo+1) * latconst > xlo) ilo--;
        ilo++;
        ihi = (xhi / latconst);
        while (ihi * latconst <= xhi) ihi++;
        ihi--;
      }

    if (is_yperiodic) {
        jlo = 0;
        jhi = ny-1;
    }
    else {
        jlo = (ylo / latconst);
        while ((jlo+1)*latconst > ylo) jlo--;
        jlo++;
        jhi = (yhi / latconst);
        while (jhi*latconst <= yhi) jhi++;
        jhi--;
    }

    if (is_zperiodic) {
        klo = 0;
        khi = nz-1;
    }
    else {
        klo = (zlo / latconst);
        while ((klo+1)*latconst > zlo) klo--;
        klo++;
        khi = (zhi / latconst);
        while (khi*latconst <= zhi) khi++;
        khi--;
    }

    
    // generate xyz coords and store them with site ID
    // tile the simulation box from origin, respecting PBC
    // site IDs should be contiguous if style = BOX and fully periodic
    // for non-periodic dims, check if site is within global box
    // for style = REGION, check if site is within region
    // if non-periodic or style = REGION, IDs may not be contiguous

    float x,y,z;
    bool can_adsorb;

    for (int k = klo; k <= khi; k++)
        for (int j = jlo; j <= jhi; j++)
            for (int i = ilo; i <= ihi; i++) {
                x = i*latconst;
                y = j*latconst;
                z = k*latconst;

                if (i == ihi && !is_xperiodic) {
                    can_adsorb = true;
                }
                else if (j == jhi && !is_yperiodic) {
                    can_adsorb = true;
                }
                else if (k == khi && !is_zperiodic) {
                    can_adsorb = true;
                }

                // By default, assume all lattice sites empty
                // TODO: This should use the global variable EMPTY_SITE
                // Don't update neighbors, since we'll use the connectivity function next
                add_site(i,j,k,x,y,z,can_adsorb,false,false);

    }


} // structered_lattice()

/* ---------------------------------------------------------------------- */

void Lattice::structured_connectivity() {

    int ineigh,jneigh,kneigh;
    uint32_t gid;
    int xneigh,yneigh,zneigh;
    
    float xprd = xhi - xlo;
    float yprd = yhi - ylo;
    float zprd = zhi - zlo;
    
    int nx = static_cast<int> (xprd / latconst);
    int ny = static_cast<int> (yprd / latconst);
    int nz = static_cast<int> (zprd / latconst);

    // create connectivity offsets
    
    int **cmap;                 // connectivity map for regular lattices
                                // cmap[maxneigh][3]
                                // 0,1,2 = i,j,k lattice unit cell offsets
    
    create(cmap,maxneigh,3,"create:cmap");

    offsets_3d(cmap);
    
    // generate global lattice connectivity for each site
    for (int i = 0; i < nsites; i++) {
        numneigh[i] = 0;

        uint32_t* neighi;
        create(neighi, maxneigh, "create:neighi");
        idneigh[i] = neighi;
      
        for (int neigh = 0; neigh < maxneigh; neigh++) {

            // ijkm neigh = indices of neighbor site
            // calculated from siteijk and cmap offsets

            ineigh = static_cast<int> (sites[i].x/latconst) + cmap[neigh][0];
            jneigh = static_cast<int> (sites[i].y/latconst) + cmap[neigh][1];
            kneigh = static_cast<int> (sites[i].z/latconst) + cmap[neigh][2];

            // xyz neigh = coords of neighbor site
            // calculated in same manner that structured_lattice() generated coords
            
            xneigh = static_cast<float> (ineigh) * latconst;
            yneigh = static_cast<float> (jneigh) * latconst;
            zneigh = static_cast<float> (kneigh) * latconst;
            
            // remap neighbor coords and indices into periodic box via ijk neigh
            // remap neighbor coords and indices into periodic box via ijk neigh

            if (is_xperiodic) {
                if (ineigh < 0) {
                    xneigh += xprd;
                    ineigh += nx;
                }
                if (ineigh >= nx) {
                    xneigh -= xprd;
                    xneigh = MAX(xneigh, static_cast<int> (xlo));
                    ineigh -= nx;
                }
            }
            if (is_yperiodic) {
                if (jneigh < 0) {
                    yneigh += yprd;
                    jneigh += ny;
                }
                if (jneigh >= ny) {
                    yneigh -= yprd;
                    yneigh = MAX(yneigh, static_cast<int> (ylo));
                    jneigh -= ny;
                }
            }
            if (is_zperiodic) {
                if (kneigh < 0) {
                    zneigh += zprd;
                    kneigh += nz;
                }
                if (kneigh >= nz) {
                    zneigh -= zprd;
                    zneigh = MAX(zneigh, static_cast<int> (zlo));
                    kneigh -= nz;
                }
            }

            // discard neighs that are outside non-periodic box or region
            if (!is_xperiodic && (xneigh < xlo || xneigh > xhi)) continue;
            if (!is_yperiodic && (yneigh < ylo || yneigh > yhi)) continue;
            if (!is_zperiodic && (zneigh < zlo || zneigh > zhi)) continue;

            // gid = global ID of neighbor
            // calculated in same manner that structured_lattice() generated IDs

            gid = ((kneigh-klo)*(jhi-jlo+1)*(ihi-ilo+1)) +
                  ((jneigh-jlo)*(ihi-ilo+1)) + ((ineigh-ilo));
            
        // add gid to neigh list of site i
        idneigh[i][numneigh[i]++] = gid;
      }
    }
    
    destroy(cmap);

} // structured_connectivity()

/* ---------------------------------------------------------------------- */

void Lattice::offsets_3d(int **cmapin) {
    
    int n = 0;
    double delx,dely,delz,r;
    double cutoff = latconst;

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {
                delx = i * latconst;
                dely = j * latconst;
                delz = k * latconst;
                r = sqrt(delx * delx + dely * dely + delz * delz);
                if (r > cutoff - EPSILON && r < cutoff + EPSILON) {
                    assert(n != maxneigh);
                    cmapin[n][0] = i;
                    cmapin[n][1] = j;
                    cmapin[n][2] = k;
                    n++;
                }
            }
        }
    }
  assert(n == maxneigh);

} // offsets_3d

/* ---------------------------------------------------------------------- */

void Lattice::add_site(uint32_t i_in, uint32_t j_in, 
                       uint32_t k_in, float x_in, float y_in, float z_in,
                       bool can_adsorb_in, bool update_neighbors_in, bool meta_neighbors_in) {
    

    std::tuple<uint32_t, uint32_t, uint32_t> key = {i_in, j_in, k_in};
    if(loc_map.find(key) != loc_map.end()) {
        // site already exists
        return;
    }

    if (nsites == nmax) {
        nmax += DELTA;
        //can_adsorb.reserve(nmax);
        //can_desorb.reserve(nmax);
    }

    // Initialize neighbor information for this new site
    numneigh[nsites] = 0;

    uint32_t* neighi;
    create(neighi, maxneigh, "create:neighi");
    idneigh[nsites]= neighi;

    // initially empty site
    sites[nsites] = Site{i_in, j_in, k_in, x_in, y_in, z_in, SPECIES_EMPTY, can_adsorb_in};
    
    loc_map[key] = nsites;
    
    if(can_adsorb_in) {
        edges[nsites] = 'a';
    }

    if (update_neighbors_in) {
        update_neighbors(nsites, meta_neighbors_in);
    }

    if (sites[nsites].x < xlo) {
        xlo = sites[nsites].x;
        ilo = sites[nsites].i;
    }
    else if (sites[nsites].x > xhi) {
        xhi = sites[nsites].x;
        ihi = sites[nsites].i;
    }

    if (sites[nsites].y < ylo) {
        ylo = sites[nsites].y;
        jlo = sites[nsites].j;
    }
    else if (sites[nsites].y > yhi) {
        yhi = sites[nsites].y;
        jhi = sites[nsites].j;
    }

    if (sites[nsites].z < zlo) {
        zlo = sites[nsites].z;
        klo = sites[nsites].k;
    }
    else if (sites[nsites].z > zhi) {
        zhi = sites[nsites].z;
        khi = sites[nsites].k;
    }

    nsites++;

    // update running max distance // TODO: make general for all types of periodicity 
    if(z_in > maxz) {
        maxz = z_in;
    }

} // add_site()

/* ---------------------------------------------------------------------- */

void Lattice::delete_site(int id) {

    loc_map.erase({sites[id].i, sites[id].j, sites[id].k});

    assert(sites.find(id) != sites.end());

    update_neighbors(id, true);

    // remove from sites 
    loc_map.erase({sites[id].i, sites[id].j, sites[id].k});
    sites.erase(id);

    // update neighbors
    for(int i = 0; i < int(numneigh[id]); i++) {
        update_neighbors(idneigh[id][i], false);
    }
    

    // delete from other hashes
    numneigh.erase(id);
    sfree(idneigh[id]);
    idneigh.erase(id);
    edges.erase(id);

} // delete_site()

/* ---------------------------------------------------------------------- */

void Lattice::update_neighbors(uint32_t n, bool meta_neighbors_in) {

    float xprd = xhi - xlo;
    float yprd = yhi - ylo;
    float zprd = zhi - zlo;
    
    uint32_t nx = static_cast<int> (xprd / latconst);
    uint32_t ny = static_cast<int> (yprd / latconst);
    uint32_t nz = static_cast<int> (zprd / latconst);

    uint32_t left, right, backward, forward, down, up;
    
    left = sites[n].i - 1;
    right = sites[n].i + 1;
    if (is_xperiodic) {
        if (left < 0) {
            left += nx;
        }
        if (right >= nx) {
            right -= nx;
        }
    }
    
    backward = sites[n].j - 1;
    forward = sites[n].j + 1;
    if (is_yperiodic) {
        if (backward < 0) {
            backward += ny;
        }
        if (forward >= ny) {
            forward -= ny;
        }
    }

    down = sites[n].k - 1;
    up = sites[n].k + 1;
    if (is_zperiodic) {
        if (down < 0) {
            down += nz;
        }
        if (up >= nz) {
            up -= nz;
        }
    }

    std::vector<std::tuple<uint32_t,uint32_t,uint32_t>> ijk;
    ijk.resize(6);

    ijk[0] = {left, sites[n].j, sites[n].k};
    ijk[1] = {right, sites[n].j, sites[n].k};
    ijk[2] = {sites[n].i, backward, sites[n].k};
    ijk[3] = {sites[n].i, forward, sites[n].k};
    ijk[4] = {sites[n].i, sites[n].j, down};
    ijk[5] = {sites[n].i, sites[n].j, up};
    

    // reset numneigh, and idneigh
    numneigh[n] = 0;
    uint32_t* neighi;
    create(neighi, maxneigh, "create:neighi");
    idneigh[n] = neighi;

    for (int q = 0; q < 6; q++) {
        auto it = loc_map.find(ijk[q]);
        if (it != loc_map.end()) {
            idneigh[n][numneigh[n]++] = it->second;
            
            if(meta_neighbors_in) {
                update_neighbors(it->second, false);
            }

        }
    }

} // update_neighbors()

/* ---------------------------------------------------------------------- */

void *Lattice::smalloc(int nbytes, const char *name)
{
    if (nbytes == 0) return NULL;

    void *ptr = malloc(size_t(nbytes));
    if(ptr == NULL) {
        std::cout << name << '\n';
        assert(false);
    }               
    return ptr;
} // smalloc()

/* ---------------------------------------------------------------------- */

template <typename TYPE>
TYPE *Lattice::create(TYPE *&array, int n, const char *name)
{
    int nbytes = ((int) sizeof(TYPE)) * n;
    array = (TYPE *) smalloc(nbytes, name);
    return array;

} // create()

/* ---------------------------------------------------------------------- */

template <typename TYPE>
TYPE **Lattice::create(TYPE **&array, int n1, int n2, const char *name) 
{
    int nbytes = ((int) sizeof(TYPE)) * n1*n2;
    TYPE *data = (TYPE *) smalloc(nbytes,name);
    nbytes = ((int) sizeof(TYPE *)) * n1;
    array = (TYPE **) smalloc(nbytes,name);
    
    int n = 0;
    for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
    }
    return array;

} // create()

/* ---------------------------------------------------------------------- */

template <typename TYPE>
void Lattice::destroy(TYPE **array)
{
    if (array == NULL) return;
    sfree(array[0]);
    sfree(array);

} // destroy()

/* ---------------------------------------------------------------------- */

void Lattice::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);

} // sfree()

/* ---------------------------------------------------------------------- */

float Lattice::get_latconst() {
    return latconst;

} // get_latconst()

/* ---------------------------------------------------------------------- */

float Lattice::get_maxz() {
    return maxz;
}

/* ---------------------------------------------------------------------- */

void Lattice::fill(std::string filename) {

    std::ifstream fin;
    fin.open(filename);

    if(!fin.is_open()) {
        std::cout << "Failed to open file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    char type;
    char junk;
    double i_in, j_in, k_in;
    int species;

    fin >> type;

    if(type == 'L') {
        while(fin >> junk >> i_in >> junk >> j_in >> junk >> k_in >> junk >> species) {
            std::tuple<uint32_t, uint32_t, uint32_t> key = {i_in, j_in, k_in};
            sites[loc_map[key]].species = species;
        }
    }
    else if(type == 'A') {
        
        for(int k = klo; k <= khi; k++) {
            fin >> junk;
            for(int i = ilo; i <= ihi; i++) {
                for(int j = jlo; j <= jhi; j++) {
                    fin >> species;
                    std::tuple<uint32_t, uint32_t, uint32_t> key = {i, j, k}; 
                    sites[loc_map[key]].species = species;
                }
            }
            fin >> junk;
        }

    }
    else {
        std::cout << "Incorrect type of input: " << type << "\n";
        exit(EXIT_FAILURE);
    }

} // fill()
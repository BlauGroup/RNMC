#include "lattice.h"
#include <cassert> 
#include "stdio.h"
#include "stdlib.h"
#include <string>
#include "math.h"
#include <iostream>
#include <map>

using namespace LGMC_NS;

#define DELTALOCAL 10000
#define DELTA 32768
#define EPSILON 0.0001


/* ---------------------------------------------------------------------- */

Lattice::Lattice(float latconst_in, 
        int ilo_in, int ihi_in, int jlo_in,
        int jhi_in, int klo_in, int khi_in, 
        bool is_xperiodic_in, bool is_yperiodic_in, bool is_zperiodic_in)  {
    // TODO: implement error handling
    
    memory = new Memory();

    latconst = latconst_in;
    
    // region of simulation input * lattice spacing
    xlo = ilo_in * latconst;
    xhi = ihi_in * latconst;
    ylo = jlo_in * latconst;
    yhi = jhi_in * latconst;
    zlo = klo_in * latconst;
    zhi = khi_in * latconst;
    
    // 0 = non-periodic, 1 = periodic
    is_xperiodic = is_xperiodic_in;
    is_yperiodic = is_yperiodic_in;
    is_zperiodic = is_zperiodic_in;
    
    nsites = 0;
    nmax = DELTA;
    sites.reserve(nmax);

    maxneigh = 6;
    idneigh.resize(0);
    numneigh.reserve(nmax);

    // create sites on lattice
    structured_lattice();
    
    // set neighbors of each site
    structured_connectivity();
    
    /*for(int n = 0; n < nsites; ++n) {
        std::cout << "id: " << n << " [" << xyz[n][0] << ", " <<
        xyz[n][1] << ", " << xyz[n][2] << "]" << std::endl;
    }*/

} // Lattice()

/* ---------------------------------------------------------------------- */

Lattice::~Lattice() {

    memory->destroy(sites);
    
    memory->destroy(numneigh);
    memory->destroy(idneigh);

    delete memory;
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

/* ----------------------------------------------------------------------
   generate site connectivity for on-lattice applications
 ------------------------------------------------------------------------- */

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
    
    memory->create(cmap,maxneigh,3,"create:cmap");

    offsets_3d(cmap);
    
    // generate global lattice connectivity for each site
    for (int i = 0; i < nsites; i++) {
        numneigh[i] = 0;

        uint32_t* neighi;
        memory->create(neighi, maxneigh, "create:neighi");
        idneigh.push_back(neighi);
      
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
            zneigh = static_cast<float> (jneigh) * latconst;
            
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
            gid = uint32_t((kneigh-klo)*(jhi-jlo+1)*(ihi-ilo+1)) +
                  uint32_t((jneigh-jlo)*(ihi-ilo+1)) + uint32_t((ineigh-ilo));
            
            
        /*std::cout << "neighbor: (" << sites[gid].x << ", " << sites[gid].y << ", " 
                  << sites[gid].z << ") for: " << "[" << sites[i].x << ", " << 
                  sites[i].y << ", " << sites[i].z << "]" << std::endl;*/
            
        // add gid to neigh list of site i
        idneigh[i][numneigh[i]++] = gid;
      }
    }

    std::cout << "numneigh" << std::endl;
    for(int i = 0; i < nsites; i++) {
        std::cout << "[" << sites[i].x << ", " <<
        sites[i].y << ", " << sites[i].z << "]" << ",";
        std::cout << "num: " << numneigh[i] << std::endl;
    }
    
    std::cout << "idneigh" << std::endl;
    for(int i = 0; i < nsites; i++) {
        std::cout << "neighbors: ";
        for(int j = 0; j < maxneigh; j++) {
            std::cout << idneigh[i][j] << ", ";
        }
        std::cout << std::endl;
    }
    
    memory->destroy(cmap);
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
    if (nsites == nmax) {
        nmax += DELTA;
        sites.reserve(nmax);
        numneigh.reserve(nmax);
        idneigh.reserve(nmax);

        // Initialize neighbor information for this new site
        numneigh.push_back(0);

        uint32_t* neighi;
        memory->create(neighi, maxneigh, "create:neighi");
        idneigh.push_back(neighi);

    }

    // initially empty site, species = 0
    sites.push_back(Site{i_in, j_in, k_in, x_in, y_in, z_in, 0, can_adsorb_in});
    std::tuple<uint32_t, uint32_t, uint32_t> key = {i_in, j_in, k_in};
    loc_map[key] = nsites;

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

} // add_site()

void Lattice::update_neighbors(uint32_t n, bool meta_neighbors_in) {

    float xprd = xhi - xlo;
    float yprd = yhi - ylo;
    float zprd = zhi - zlo;
    
    int nx = static_cast<int> (xprd / latconst);
    int ny = static_cast<int> (yprd / latconst);
    int nz = static_cast<int> (zprd / latconst);

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
            backward += nx;
        }
        if (forward >= nx) {
            forward -= nx;
        }
    }

    down = sites[n].k - 1;
    up = sites[n].k + 1;
    if (is_zperiodic) {
        if (down < 0) {
            down += nx;
        }
        if (up >= nx) {
            up -= nx;
        }
    }

    std::tuple<uint32_t,uint32_t,uint32_t> *ijk = {
        {left, sites[n].j, sites[n].k},
        {right, sites[n].j, sites[n].k},
        {sites[n].i, backward, sites[n].k},
        {sites[n].i, forward, sites[n].k},
        {sites[n].i, sites[n].j, down},
        {sites[n].i, sites[n].j, up},
    };

    std::map<std::tuple<uint32_t,uint32_t,uint32_t>, int>::iterator it;

    int thisnumneigh = 0;
    for (int q = 0; q < 6; q++) {
        it = loc_map.find(ijk[q]);
        if (it != loc_map.end()) {
            idneigh[n][thisnumneigh] = it->second;

            if(meta_neighbors_in) {
                update_neighbors(idneigh[n][thisnumneigh], false);
            }

            thisnumneigh++;
        }
    }

} // update_neighbors()
#include "lattice.h"
#include <cassert> 
#include "stdio.h"
#include "stdlib.h"
#include "string"
#include "math.h"
#include <numeric>
#include "../GMC/GMC.cpp"
#include <iostream>
#include <format>

using namespace LGMC_NS;

#define DELTALOCAL 10000
#define DELTA 32768
#define EPSILON 0.0001
#define DELTABUF 10000

// types of reactions
enum {ELECTRO, DIFFU, CHEM, DESORP};


/* ---------------------------------------------------------------------- */

Lattice::Lattice(int latconst_in, <Solver, Model, Parameters, TrajectoriesSql>Dispatcher *ptr_in, 
        int boxxlo_in, int boxxhi_in, int boxylo_in,
        int boxyhi_in, int boxzlo_in, int boxzhi_in)  {
    // TODO: implement error handling
    
    memory = new Memory();

    latconst = latconst_in;

    // pointer to gillespie dispatcher
    dis_ptr = ptr_in;
    
    // region of simulation input * lattice spacing
    boxxlo = boxxlo_in;
    boxxhi = boxxhi_in;
    boxylo = boxylo_in;
    boxyhi = boxyhi_in;
    boxzlo = boxzlo_in;
    boxzhi = boxzhi_in;
    
    // 0 = non-periodic, 1 = periodic
    xperiodic = 1;
    yperiodic = 1;
    zperiodic = 0;
    
    nsites = nmax = 0;
    sites = NULL;
    
    maxneigh = 6;
    numneigh = NULL;
    idneigh = NULL;
    
    prop_sum = 0;
    
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
    
    // set domain->nx,ny,nz iff style = BOX and system is fully periodic
    // else site IDs may be non-contiguous and/or ordered irregularly
    // 3 dimensions
    uint32_t nx, ny, nz;
    nx = (boxxhi - boxxlo / latconst);
    ny = (boxyhi - boxylo / latconst);
    nz = (boxzhi - boxzlo / latconst);

    // if dim is periodic:
    //    lattice origin = lower box boundary
    //    loop bounds = 0 to N-1
    // if dim is non-periodic:
    //   lattice origin = 0.0
    //   loop bounds = enough to tile box completely, with all basis atoms
    
    if (xperiodic) {
        xlo = 0;
        xhi = nx-1;
      }
    else {
        xlo = (boxxlo / latconst);
        while ((xlo+1)*latconst > boxxlo) xlo--;
        xlo++;
        xhi = (boxxhi / latconst);
        while (xhi*latconst <= boxxhi) xhi++;
        xhi--;
      }

    if (yperiodic) {
        ylo = 0;
        yhi = ny-1;
    }
    else {
        ylo = (boxylo / latconst);
        while ((ylo+1)*latconst > boxylo) ylo--;
        ylo++;
        yhi = (boxyhi / latconst);
        while (yhi*latconst <= boxyhi) yhi++;
        yhi--;
    }

    if (zperiodic) {
        zlo = 0;
        zhi = nz-1;
    }
    else {
        zlo = (boxzlo / latconst);
        while ((zlo+1)*latconst > boxzlo) zlo--;
        zlo++;
        zhi = (boxzhi / latconst);
        while (zhi*latconst <= boxzhi) zhi++;
        zhi--;
    }

    
    // generate xyz coords and store them with site ID
    // tile the simulation box from origin, respecting PBC
    // site IDs should be contiguous if style = BOX and fully periodic
    // for non-periodic dims, check if site is within global box
    // for style = REGION, check if site is within region
    // if non-periodic or style = REGION, IDs may not be contiguous

    uint32_t x,y,z;

    uint32_t n = 0;
    for (int k = zlo; k <= zhi; k++)
        for (int j = ylo; j <= yhi; j++)
            for (int i = xlo; i <= xhi; i++) {
                n++;
                x = i*latconst;
                y = j*latconst;
                z = k*latconst;

                add_site(n,x,y,z);

    }


} // structered_lattice()

/* ----------------------------------------------------------------------
   generate site connectivity for on-lattice applications
 ------------------------------------------------------------------------- */

void Lattice::structured_connectivity() {

    int ineigh,jneigh,kneigh;
    int gid;
    int xneigh,yneigh,zneigh;
    int xprd, yprd, zprd;
    
    xprd = boxxhi - boxxlo;
    yprd = boxyhi - boxylo;
    zprd = boxzhi - boxzlo;
    
    int nx, ny, nz;
    nx = xprd / latconst;
    ny = yprd / latconst;
    nz = zprd / latconst;
    
    memory->create(idneigh,nsites,maxneigh,"create:idneigh");
    memory->create(numneigh,nmax,"create:numneigh");
    
    // create connectivity offsets
    
    int **cmap;                 // connectivity map for regular lattices
                                // cmap[maxneigh][3]
                                // 0,1,2 = i,j,k lattice unit cell offsets
    
    memory->create(cmap,maxneigh,3,"create:cmap");

    offsets_3d(cmap);
    
    // generate global lattice connectivity for each site
    for (int i = 0; i < nsites; i++) {
        numneigh[i] = 0;
      
        for (int neigh = 0; neigh < maxneigh; neigh++) {

            // ijkm neigh = indices of neighbor site
            // calculated from siteijk and cmap offsets

            ineigh = static_cast<int> (sites[i].x/latconst) + cmap[neigh][0];
            jneigh = static_cast<int> (sites[i].y/latconst) + cmap[neigh][1];
            kneigh = static_cast<int> (sites[i].z/latconst) + cmap[neigh][2];

            // xyz neigh = coords of neighbor site
            // calculated in same manner that structured_lattice() generated coords
            
            xneigh = ineigh * static_cast<int> (latconst);
            yneigh = jneigh * static_cast<int> (latconst);
            zneigh = kneigh * static_cast<int> (latconst);
            
            // remap neighbor coords and indices into periodic box via ijk neigh
            // remap neighbor coords and indices into periodic box via ijk neigh

            if (xperiodic) {
                if (ineigh < 0) {
                    xneigh += xprd;
                    ineigh += nx;
                }
                if (ineigh >= nx) {
                    xneigh -= xprd;
                    xneigh = MAX(xneigh, static_cast<int> (boxxlo));
                    ineigh -= nx;
                }
            }
            if (yperiodic) {
                if (jneigh < 0) {
                    yneigh += yprd;
                    jneigh += ny;
                }
                if (jneigh >= ny) {
                    yneigh -= yprd;
                    yneigh = MAX(yneigh, static_cast<int> (boxylo));
                    jneigh -= ny;
                }
            }
            if (zperiodic) {
                if (kneigh < 0) {
                    zneigh += zprd;
                    kneigh += nz;
                }
                if (kneigh >= nz) {
                    zneigh -= zprd;
                    zneigh = MAX(zneigh, static_cast<int> (boxzlo));
                    kneigh -= nz;
                }
            }

            // discard neighs that are outside non-periodic box or region
            if (!xperiodic && (xneigh < boxxlo || xneigh > boxxhi)) continue;
            if (!yperiodic && (yneigh < boxylo || yneigh > boxyhi)) continue;
            if (!zperiodic && (zneigh < boxzlo || zneigh > boxzhi)) continue;

            // gid = global ID of neighbor
            // calculated in same manner that structured_lattice() generated IDs
            uint32_t one = 1;   // use this to avoid int overflow in calc of gid
            gid = one*uint32_t((kneigh-zlo)*(yhi-ylo+1)*(xhi-xlo+1)) +
                  one*uint32_t((jneigh-ylo)*(xhi-xlo+1)) + one*uint32_t((ineigh-xlo));
            
            
        std::cout << "neighbor: (" << sites[gid].x << ", " << sites[gid].y << ", " 
                  << sites[gid].z << ") for: " << "[" << sites[i].x << ", " << 
                  sites[i].y << ", " << sites[i].z << "]" << std::endl;
            
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
    double cutoff = 1*latconst;

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

void Lattice::add_site(uint32_t n, uint32_t x, uint32_t y, uint32_t z) {
    if (nsites == nmax) grow(0);

    // initially empty site, species = -1
    sites[n] = Site{x, y, z, -1};

    nsites++;
} // add_site()

/* ---------------------------------------------------------------------- */

void Lattice::grow(uint32_t n) {
    if (n == 0) nmax += DELTA;
    else nmax = n;

    memory->grow(sites,nmax,"grow:sites");
    
} // grow()

/* ---------------------------------------------------------------------- 
    Only calls this function if necessary reactants are on lattice sites
---------------------------------------------------------------------- */

void Lattice::update_propensity(int site_one, int site_two, Lat_Reaction &reaction) {
    
    // TODO: for certain reactions make sure there is an empty site 
    // TODO: for interactions with electrolyte make sure possible (take into account distance??)

    double p;
    // zero reactants
    if (reaction.number_of_reactants == 0)
        p = dis_ptr->model.factor_zero * reaction.rate;

    // one reactant
    else if (reaction.number_of_reactants == 1)
        p = 1 * reaction.rate;


    // two reactants
    else {
       if (reaction.reactants[0] == reaction.reactants[1])
            p = dis_ptr->model.factor_duplicate
                * dis_ptr->model.factor_two
                * 2
                * (2 - 1)
                * reaction.rate;

        else
            p = dis_ptr->model.factor_two
                * 1
                * 1
                * reaction.rate;
    }

    // add or change existing propensity 
    std::string site_combo = (site_one < site_two) ? std::to_string(site_one) + "." + std::to_string(site_two) : std::to_string(site_two) + "." + std::to_string(site_one);
    props[site_combo].push_back(p);

    // update running sum 
    prop_sum += p;


    
} // update_propensity()

/* ---------------------------------------------------------------------- */

void Lattice::update(int site_one, int site_two) {

    if(site_one > -1) {
        relevant_react(site_one);       
    }
    if(site_two > -1) {
        relevant_react(site_two);
    }
    
}

/* ---------------------------------------------------------------------- */

void Lattice::clear_site(int site) {

    // reset or initiate site combos
    for(uint32_t neigh = 0; neigh < numneigh[site]; neigh++ ) {
        int neighbor = idneigh[site][neigh];

        std::string combo = (site > neighbor) ? std::to_string(site) + "." + std::to_string(neighbor) : std::to_string(neighbor) + "." + std::to_string(site);
        
        // check if first time key has been added
        if(props.find(combo) == props.end()) {
            // key not found
            // TODO: is this valid?????
            props[combo].resize(lat_dependents.size());
        }
        else {
            // already exists, clear vector to update
            prop_sum -= std::accumulate(props[combo].begin(), props[combo].end(), 0);

            // TODO: make sure capacity does not get changed
            props[combo].clear();
        }
    }

    // reset for gillespie (-2) and empty site (-1)
    for(int i = -2; i < 0; i++) {
        std::string combo = std::to_string(site) + "." + std::to_string(i);

        if(props.find(combo) == props.end()) {
         props[combo].resize(lat_dependents.size());
        }
        else {
            prop_sum -= std::accumulate(props[combo].begin(), props[combo].end(), 0);
            props[combo].clear();
        }

    }

} // clear_site

/* ---------------------------------------------------------------------- */

void Lattice::relevant_react(int site) {

    clear_site(site);

    // all reactions related to central site 
    std::vector<int> &potential_reactions = lat_dependents[sites[site].species]; 

    // compute and add new propensities 
    for(int reaction_id = 0; reaction_id < static_cast<int> (potential_reactions.size()); reaction_id++ ) {

        Lat_Reaction &reaction = dis_ptr->model.reactions[reaction_id];

        // reaction with gillepsie 
        if(reaction.type == DESORP) {
            update_propensity(site, -2, reaction);
        }

        // single reaction 
        if(reaction.number_of_reactants == 1) {
            
            if(reaction.number_of_products == 1) {
                update_propensity(site, -1, reaction);
            }
            else {
                // two products make sure site is empty
                for(uint32_t neigh = 0; neigh < numneigh[site]; neigh++) {
                    
                    int neighbor = idneigh[site][neigh];
                    
                    if(sites[neighbor].species == -1) {
                        // empty 
                        update_propensity(site, neighbor, reaction);
                    }


                } // for neigh
            }

        } // single reactant

        if(reaction.number_of_reactants == 2) {
    
            int other_reactant = (site == reaction.reactants[0]) ? reaction.reactants[1] : reaction.reactants[0];
            // make sure neighbor is relevant 
            for(uint32_t neigh = 0; neigh < numneigh[site]; neigh++) {
                int neighbor = idneigh[site][neigh];
                    
                if(sites[neighbor].species == other_reactant) {
                    update_propensity(site, neighbor, reaction);
                }

            } // for neigh
            
        } // two reactants


    }

}

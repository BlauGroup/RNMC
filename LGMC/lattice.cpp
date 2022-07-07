#include "lattice.h"
#include "stdio.h"
#include "stdlib.h"
#include "string"
#include "math.h"
#include <iostream>

using namespace LGMC_NS;

#define DELTALOCAL 10000
#define DELTA 32768
#define EPSILON 0.0001
#define DELTABUF 10000

/* ---------------------------------------------------------------------- */

Lattice::Lattice(int, char **arg) {
    // TODO: implement error handling
    
    memory = new Memory();

    latconst = atof(arg[1]);
    
    // region of simulation input * lattice spacing
    boxxlo = latconst*atof(arg[2]);
    boxxhi = latconst*atof(arg[3]);
    boxylo = latconst*atof(arg[4]);
    boxyhi = latconst*atof(arg[5]);
    boxzlo = latconst*atof(arg[6]);
    boxzhi = latconst*atof(arg[7]);
    
    // 0 = non-periodic, 1 = periodic
    xperiodic = 1;
    yperiodic = 1;
    zperiodic = 0;
    
    nsites = nmax = 0;
    id = NULL;
    xyz = NULL;
    
    maxneigh = 6;
    numneigh = NULL;
    idneigh = NULL;
    
    // create sites on lattice
    structured_lattice();
    
    // set neighbors of each site
    structured_connectivity();

    // TODO: set initial state and propensities
    
    /*for(int n = 0; n < nsites; ++n) {
        std::cout << "id: " << n << " [" << xyz[n][0] << ", " <<
        xyz[n][1] << ", " << xyz[n][2] << "]" << std::endl;
    }*/

} // Lattice()

/* ---------------------------------------------------------------------- */

Lattice::~Lattice() {

    memory->destroy(id);
    memory->destroy(xyz);
    
    memory->destroy(numneigh);
    memory->destroy(idneigh);

    delete memory;
} // ~Lattice()


/* ---------------------------------------------------------------------- */

void Lattice::structured_lattice() {
    
    // set domain->nx,ny,nz iff style = BOX and system is fully periodic
    // else site IDs may be non-contiguous and/or ordered irregularly
    // 3 dimensions
    int nx, ny, nz;
    nx = static_cast<int> (boxxhi - boxxlo / latconst);
    ny = static_cast<int> (boxyhi - boxylo / latconst);
    nz = static_cast<int> (boxzhi - boxzlo / latconst);

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
        xlo = static_cast<int> (boxxlo / latconst);
        while ((xlo+1)*latconst > boxxlo) xlo--;
        xlo++;
        xhi = static_cast<int> (boxxhi / latconst);
        while (xhi*latconst <= boxxhi) xhi++;
        xhi--;
      }

    if (yperiodic) {
        ylo = 0;
        yhi = ny-1;
    }
    else {
        ylo = static_cast<int> (boxylo / latconst);
        while ((ylo+1)*latconst > boxylo) ylo--;
        ylo++;
        yhi = static_cast<int> (boxyhi / latconst);
        while (yhi*latconst <= boxyhi) yhi++;
        yhi--;
    }

    if (zperiodic) {
        zlo = 0;
        zhi = nz-1;
    }
    else {
        zlo = static_cast<int> (boxzlo / latconst);
        while ((zlo+1)*latconst > boxzlo) zlo--;
        zlo++;
        zhi = static_cast<int> (boxzhi / latconst);
        while (zhi*latconst <= boxzhi) zhi++;
        zhi--;
    }

    
    // generate xyz coords and store them with site ID
    // tile the simulation box from origin, respecting PBC
    // site IDs should be contiguous if style = BOX and fully periodic
    // for non-periodic dims, check if site is within global box
    // for style = REGION, check if site is within region
    // if non-periodic or style = REGION, IDs may not be contiguous

    double x,y,z;

    tagint n = 0;
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
    tagint gid;
    double xneigh,yneigh,zneigh;
    double xprd, yprd, zprd;
    
    xprd = boxxhi - boxxlo;
    yprd = boxyhi - boxylo;
    zprd = boxzhi - boxzlo;
    
    int nx, ny, nz;
    nx = static_cast<int> (xprd / latconst);
    ny = static_cast<int> (yprd / latconst);
    nz = static_cast<int> (zprd / latconst);
    
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

            ineigh = xyz[i][0]/latconst + cmap[neigh][0];
            jneigh = xyz[i][1]/latconst + cmap[neigh][1];
            kneigh = xyz[i][2]/latconst + cmap[neigh][2];

            // xyz neigh = coords of neighbor site
            // calculated in same manner that structured_lattice() generated coords
            
            xneigh = ineigh * latconst;
            yneigh = jneigh * latconst;
            zneigh = kneigh * latconst;
            
            // remap neighbor coords and indices into periodic box via ijk neigh
            // remap neighbor coords and indices into periodic box via ijk neigh

            if (xperiodic) {
                if (ineigh < 0) {
                    xneigh += xprd;
                    ineigh += nx;
                }
                if (ineigh >= nx) {
                    xneigh -= xprd;
                    xneigh = MAX(xneigh,boxxlo);
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
                    yneigh = MAX(yneigh,boxylo);
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
                    zneigh = MAX(zneigh,boxzlo);
                    kneigh -= nz;
                }
            }

            // discard neighs that are outside non-periodic box or region
            if (!xperiodic && (xneigh < boxxlo || xneigh > boxxhi)) continue;
            if (!yperiodic && (yneigh < boxylo || yneigh > boxyhi)) continue;
            if (!zperiodic && (zneigh < boxzlo || zneigh > boxzhi)) continue;

            // gid = global ID of neighbor
            // calculated in same manner that structured_lattice() generated IDs
            int one = 1;   // use this to avoid int overflow in calc of gid
            gid = one*tagint((kneigh-zlo)*(yhi-ylo+1)*(xhi-xlo+1)) +
                  one*tagint((jneigh-ylo)*(xhi-xlo+1)) + one*tagint((ineigh-xlo));
            
            
        std::cout << "neighbor: (" << xyz[gid][0] << ", " << xyz[gid][1] << ", " << xyz[gid][2] << ") for: " << "[" << xyz[i][0] << ", " << xyz[i][1] << ", " << xyz[i][2] << "]" << std::endl;
            
        // add gid to neigh list of site i
        idneigh[i][numneigh[i]++] = gid;
      }
    }

    // delete siteijk and connectivity offsets
    //memory->destroy(cmap);
    //memory->destroy(basis);
    std::cout << "numneigh" << std::endl;
    for(int i = 0; i < nsites; i++) {
        std::cout << "[" << xyz[i][0] << ", " <<
        xyz[i][1] << ", " << xyz[i][2] << "]" << ",";
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

void Lattice::add_site(tagint n, double x, double y, double z) {
    if (nsites == nmax) grow(0);

    id[nsites] = n;
    xyz[nsites][0] = x;
    xyz[nsites][1] = y;
    xyz[nsites][2] = z;

    nsites++;
} // add_site()

/* ---------------------------------------------------------------------- */

void Lattice::grow(int n) {
    if (n == 0) nmax += DELTA;
    else nmax = n;

    memory->grow(id,nmax,"grow:id");
    memory->grow(xyz,nmax,3,"grow:xyz");
    
} // grow()

/* ----------------------------------------------------------------------
Set species list to -1, state list to 0, set propensities to zero
 ---------------------------------------------------------------------- */

void Lattice::init() {
    
    //memory->create(state,,"create:cmap");
    
    
} // init()

/* ---------------------------------------------------------------------- */

void Lattice::update() {
    
    //memory->create(state,,"create:cmap");
    
    
} // update()

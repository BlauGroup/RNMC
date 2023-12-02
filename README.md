<img src="./logo.png">

Reaction Network Monte Carlo (`RNMC`) is a collection of programs for kinetic Monte Carlo (kMC) simulation of physical systems heavily based on [SPPARKS](https://spparks.sandia.gov/). `RNMC` is designed to run large numbers of simulations of a fixed system in parallel. The project consists of three kMC modules for different domains which use `core` code for shared processes, for example IO, threading logic and model independent simulation logic.

## Three Modules
- `GMC` - [Gillespie Monte Carlo](https://lzichi.github.io/RNMC/GMC.html): Implementation of Gillespie's next reaction simulator. `GMC` is able to run simulations of reaction networks with hundreds of millions of reactions.
- `NPMC` - [Nano Particle Monte Carlo](https://lzichi.github.io/RNMC/NPMC.html): A three dimensional statistical field theory simulator which supports one- and two-site interactions. Useful for simulating nano particles.
- `LGMC` - [Lattice Gillespie Monte Carlo](https://lzichi.github.io/RNMC/LGMC.html): A kMC implementation coupling a homogeneous (Gillespie-like) region with a lattice, enabling simulations with reactions occurring in multiple phases and capable of electrochemical reactions.

See [this](https://doi.org/10.26434/chemrxiv-2021-c2gp3) paper for an example of the kind of work being done with RNMC.

Please see the [documentation](https://lzichi.github.io/RNMC/) for more information.

## Dependencies

RNMC depends on [GSL](https://www.gnu.org/software/gsl/) for pseudo random number generation and [sqlite](https://www.sqlite.org/index.html) for the database interfaces.

## Building

On a machine with system versions of GSL and sqlite, the executables can be built with a makefile. There are makefiles inside the `GMC`, `NPMC`, or `LGMC` folders.

To make an executable for `GMC`, `NPMC` or `LGMC` first enter that folder to use the makefile. To create an executable, use make and then the name of the module. 

For example, to use `GMC`:

```
$ cd GMC
```

```
$ make GMC
```

For further help on the makefile and to view other commands:

```
$ make help
```

Note that the makefile uses the `gsl-config` utility to find headers and libraries for GSL. If you are on a cluster and sqlite is not present, it can be built as follows:

```
cd $HOME
wget https://www.sqlite.org/2021/sqlite-amalgamation-3360000.zip
unzip sqlite-amalgamation-3360000.zip
cd sqlite-amalgamation-3360000
gcc -o libsqlite3.so -shared -fPIC sqlite3.c -lpthread -ldl
```

in which case, the simulators can be built like this:

```
export CPATH=$HOME/sqlite-amalgamation-3360000:$CPATH
export LIBRARY_PATH=$HOME/sqlite-amalgamation-3360000:$LIBRARY_PATH
make GMC
```

If you need to build [GSL](https://www.gnu.org/software/gsl/) from source: 

```
wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
mkdir gsl
mv gsl-latest.tar.gz gsl
tar -xvf gsl-latest.tar.gz
cd gsl-2.7.1
./configure --prefix=/my/path/to/gsl
make
make install
echo $PKG_CONFIG_PATH 
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/my/path/to/gsl-2.7.1
```

Note that if you build from source use `pkg-config gsl` instead of `gsl-config` inside each makefile


## Testing

Inside the `tests` folder are unit tests using [GoogleTest](https://google.github.io/googletest/primer.html) and end-to-end tests for GMC and NPMC (LGMC trajectories are not deterministic). Unit tests can be run with the makefile and end-to-end tests can be run using `test.sh`.

## Running GMC/NPMC/LGMC

## GMC

GMC is run as follows:

```
$ cd GMC
$ make GMC
```

When running `GMC` ensure that your input file paths are correct considering the executable is inside the `GMC` folder. Below is an example of how `GMC` can be run using the input files inside the <a href="{{ site.github.repository_url }}"> examples directory </a> (here `step_cutoff` is specified):

```
./GMC --reaction_database=../examples/GMC/rn.sqlite --initial_state_database=../examples/GMC/initial_state.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200
```

`GMC` requires six input arguments (either `step_cutoff` or `time_cutoff` must be specified): 

- <span style="color:#0066CC"> reaction_database </span>: a sqlite database containing the reaction network and metadata.
- <span style="color:#0066CC"> initial_state_database </span>: a sqlite database containing initial state. The simulation trajectories are also written into the database.
-  <span style="color:#0066CC">number_of_simulation </span>: an integer specifying how many simulations to run.
-  <span style="color:#0066CC">base_seed </span>: seeds used are `base_seed, base_seed+1, ..., base_seed+number_of_simulations-1`.
- <span style="color:#0066CC"> thread_count </span>: how many threads to use.
- <span style="color:#0066CC"> step_cutoff </span>: how many steps in each simulation.
- <span style="color:#0066CC"> time_cutoff </span>: how much time in each simulation [s].

### The Reaction Network Database 
There are two tables in the reaction network database both of which **must be created and filled in by the user**:

```
CREATE TABLE metadata (
        number_of_species   INTEGER NOT NULL,
        number_of_reactions INTEGER NOT NULL
);
```

```
CREATE TABLE reactions (
        reaction_id         INTEGER NOT NULL PRIMARY KEY,
        number_of_reactants INTEGER NOT NULL,
        number_of_products  INTEGER NOT NULL,
        reactant_1          INTEGER NOT NULL,
        reactant_2          INTEGER NOT NULL,
        product_1           INTEGER NOT NULL,
        product_2           INTEGER NOT NULL,
        rate                REAL NOT NULL
);
```

### The State Database 
There are five tables in the initial state database all of which **must be created by the user**: 

```
CREATE TABLE initial_state (
        species_id             INTEGER NOT NULL PRIMARY KEY,
        count                  INTEGER NOT NULL
);
```
```
CREATE TABLE trajectories (
        seed                INTEGER NOT NULL,
        step                INTEGER NOT NULL,
        time                REAL NOT NULL,
        reaction_id         INTEGER NOT NULL
);
```
```
CREATE TABLE factors (
        factor_zero      REAL NOT NULL,
        factor_two       REAL NOT NULL,
        factor_duplicate REAL NOT NULL
);
```
```
CREATE TABLE interrupt_state (
        seed                    INTEGER NOT NULL,
        species_id              INTEGER NOT NULL,
        count                   INTEGER NOT NULL
);
```
```
CREATE TABLE interrupt_cutoff (
        seed                    INTEGER NOT NULL,
        step                    INTEGER NOT NULL,
        time                    INTEGER NOT NULL       
);
```

## Running NPMC
NPMC is run as follows:

```
$ cd NPMC
$ make NPMC
```

```
./NPMC --nano_particle_database=../examples/NPMC/np.sqlite --initial_state_database=../examples/NPMC/initial_state.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200 
```

`NPMC` requires six input arguments (either `step_cutoff` or `time_cutoff` must be specified): 

- <span style="color:#0066CC"> nano_particle_database </span>: a sqlite database containing the nanoparticle data and metadata.
- <span style="color:#0066CC"> initial_state_database </span> : a sqlite database containing initial state. The simulation trajectories are also written into the database
- <span style="color:#0066CC"> number_of_simulation </span>: an integer specifying how many simulations to run
-  <span style="color:#0066CC"> base_seed </span>: seeds used are `base_seed, base_seed+1, ..., base_seed+number_of_simulations-1`
- <span style="color:#0066CC"> thread_count </span>: how many threads to use.
- <span style="color:#0066CC"> step_cutoff </span>: how many steps in each simulation.
- <span style="color:#0066CC"> time_cutoff </span>: how much time in each simulation [s].

### The Nano Particle Database
There are four tables in the nanoparticle database all of which **must be created and filled in by the user**:

```
CREATE TABLE species (
    species_id          INTEGER NOT NULL PRIMARY KEY,
    degrees_of_freedom  INTEGER NOT NULL
);
```

```
CREATE TABLE sites (
    site_id             INTEGER NOT NULL PRIMARY KEY,
    x                   REAL NOT NULL,
    y                   REAL NOT NULL,
    z                   REAL NOT NULL,
    species_id          INTEGER NOT NULL
);
```

```
CREATE TABLE interactions (
    interaction_id      INTEGER NOT NULL PRIMARY KEY,
    number_of_sites     INTEGER NOT NULL,
    species_id_1        INTEGER NOT NULL,
    species_id_2        INTEGER NOT NULL,
    left_state_1        INTEGER NOT NULL,
    left_state_2        INTEGER NOT NULL,
    right_state_1       INTEGER NOT NULL,
    right_state_2       INTEGER NOT NULL,
    rate                REAL NOT NULL
);
```

```
CREATE TABLE metadata (
    number_of_species                   INTEGER NOT NULL,
    number_of_sites                     INTEGER NOT NULL,
    number_of_interactions              INTEGER NOT NULL
);
```

### The State Database
There are five tables in the initial state database all of which **must be created by the user**: 

```
CREATE TABLE initial_state (
    site_id            INTEGER NOT NULL PRIMARY KEY,
    degree_of_freedom  INTEGER NOT NULL
);
```

```
CREATE TABLE trajectories (
    seed               INTEGER NOT NULL,
    step               INTEGER NOT NULL,
    time               REAL NOT NULL,
    site_id_1          INTEGER NOT NULL,
    site_id_2          INTEGER NOT NULL,
    interaction_id     INTEGER NOT NULL
);
```

```
CREATE TABLE factors (
    one_site_interaction_factor      REAL NOT NULL,
    two_site_interaction_factor      REAL NOT NULL,
    interaction_radius_bound         REAL NOT NULL,
    distance_factor_type             TEXT NOT NULL
);
```

```
CREATE TABLE interrupt_state (
    seed                    INTEGER NOT NULL,
    site_id                 INTEGER NOT NULL,
    degree_of_freedom       INTEGER NOT NULL
); 
```

```
CREATE TABLE interrupt_cutoff (
        seed                    INTEGER NOT NULL,
        step                    INTEGER NOT NULL,
        time                    INTEGER NOT NULL   
);
```

## LGMC
LGMC is run as follows:

```
$ cd LGMC
$ make LGMC
```

```
./LGMC --lattice_reaction_database=../examples/LGMC/CO_oxidation/rn.sqlite --initial_state_database=../examples/LGMC/CO_oxidation/initial_state.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200 --parameters=../examples/LGMC/CO_oxidation/LGMC_params.txt
```

LGMC requires seven input arguments (either `step_cutoff` or `time_cutoff` must be specified): 

- <span style="color:#0066CC"> reaction_database </span>: a sqlite database containing the reaction network and metadata.
- <span style="color:#0066CC"> initial_state_database </span>: a sqlite database containing initial state. The simulation trajectories are also written into the database.
-  <span style="color:#0066CC">number_of_simulation </span>: an integer specifying how many simulations to run.
-  <span style="color:#0066CC">base_seed </span>: seeds used are `base_seed, base_seed+1, ..., base_seed+number_of_simulations-1`.
- <span style="color:#0066CC"> thread_count </span>: how many threads to use.
- <span style="color:#0066CC"> step_cutoff </span>: how many steps in each simulation.
- <span style="color:#0066CC"> time_cutoff </span>: how much time in each simulation [s].
- <span style="color:#0066CC"> parameters </span>: name of .txt file containing additional parameters
    - <span style="color:#006633"> lattice_constant </span>
    - <span style="color:#006633"> box x upper boundary </span>
    - <span style="color:#006633"> box y upper boundary </span>
    - <span style="color:#006633"> box z upper boundary </span>
    - <span style="color:#006633"> temperature </span>
    - <span style="color:#006633"> Electron free energy </span>
    - <span style="color:#006633"> Is add site (T\|F) </span>
    - <span style="color:#006633"> Charge transfer style (M\|B) </span>

### The Lattice Reaction Network Database 

There are two tables in the lattice reaction network database both of which **must be created and filled in by the user**:

```
CREATE TABLE metadata (
    number_of_species   INTEGER NOT NULL,
    number_of_reactions INTEGER NOT NULL
);
```

```
CREATE TABLE reactions (
        reaction_id                     INTEGER NOT NULL PRIMARY KEY,
        number_of_reactants             INTEGER NOT NULL,
        number_of_products              INTEGER NOT NULL,
        reactant_1                      INTEGER NOT NULL,
        reactant_2                      INTEGER NOT NULL,
        product_1                       INTEGER NOT NULL,
        product_2                       INTEGER NOT NULL,
        phase_reactant_1                CHAR(1) NOT NULL,
        phase_reactant_2                CHAR(1) NOT NULL,
        phase_product_1                 CHAR(1) NOT NULL,
        phase_product_2                 CHAR(1) NOT NULL,
        dG                              REAL NOT NULL,
        prefactor                       REAL NOT NULL,
        rate                            REAL NOT NULL,
        electron_tunneling_coefficient  REAL NOT NULL,
        reorganization_energy           REAL NOT NULL,
        charge_transfer_coefficient     REAL NOT NULL,
        type                            CHAR(1) NOT NULL
);
```
### The State Database 
There are five tables in the initial state database all of which **must be created by the user**: 
```
CREATE TABLE initial_state (
        species_id             INTEGER NOT NULL PRIMARY KEY,
        count                  INTEGER NOT NULL
);
```
```
CREATE TABLE trajectories (
        seed                INTEGER NOT NULL,
        step                INTEGER NOT NULL,
        time                REAL NOT NULL,
        reaction_id         INTEGER NOT NULL,
        site_1_mapping      INTEGER NOT NULL,
        site_2_mapping      INTEGER NOT NULL
);
```
```
CREATE TABLE factors (
        factor_zero      REAL NOT NULL,
        factor_two       REAL NOT NULL,
        factor_duplicate REAL NOT NULL
);
```
```
CREATE TABLE interrupt_state (
        seed                    INTEGER NOT NULL,
        species_id              INTEGER NOT NULL,
        quantity                INTEGER NOT NULL,
        site_mapping            INTEGER NOT NULL,
        edge                    INTEGER NOT NULL
        
);
```
```
CREATE TABLE interrupt_cutoff (
        seed                    INTEGER NOT NULL,
        step                    INTEGER NOT NULL,
        time                    INTEGER NOT NULL,
        maxk                    INTEGER NOT NULL
);
```
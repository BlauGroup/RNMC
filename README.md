<img src="./logo.png">

Reaction Network Monte Carlo (RNMC) is a collection of programs for Monte Carlo simulation of statistical mechanical systems heavily inspired by [SPPARKS](https://spparks.sandia.gov/). RNMC is designed to run large numbers of simulations of a fixed system in parallel. The project currently consists of three parts:
- `core` : Core code shared by all simulators, for example IO, threading logic and model independent simulation logic.
- `GMC` : Implementation of Gillespie's next reaction simulator. GMC is able to run simulations of reaction networks with hundreds of millions of reactions, even when the number of species is small.
- `NPMC` : A 3D statistical field theory simulator which supports one and two site interactions. Useful for simulating nano particles.

See [this](https://doi.org/10.26434/chemrxiv-2021-c2gp3) paper for an example of the kind of work being done with RNMC.

### Dependencies

RNMC depends on [GSL](https://www.gnu.org/software/gsl/) for pseudo random number generation and [sqlite](https://www.sqlite.org/index.html) for the database interfaces.

### Building

On a machine with system versions of GSL and sqlite, the executables can be built like this:
```
CC=g++ ./build.sh
```
The executables are put in the `build` directory. Note that the build script uses the `gsl-config` utility to find headers and libraries for GSL. If you are on a cluster and sqlite is not present, it can be built as follows:

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
CC=g++ ./build.sh
```

### Testing

Run the tests using `test.sh` from the root directory of the repository.

## Running GMC

GMC is run as follows:

```
GMC --reaction_database=rn.sqlite --initial_state_database=initial_state.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200
```

- `reaction_database`: a sqlite database containing the reaction network and metadata.
- `initial_state_database` : a sqlite database containing initial state. The simulation trajectories are also written into the database
- `number_of_simulation`: an integer specifying how many simulations to run
- `base_seed`: seeds used are `base_seed, base_seed+1, ..., base_seed+number_of_simulations-1`
- `thread_count`: is how many threads to use.
- `step_cutoff`: how many steps in each simulation

### The Reaction Network Database

There are 2 tables in the reaction network database:
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
There are 3 tables in the initial state database. The factors can be used to modify rates of reactions which have zero or two reactants, or have duplicate reactants.
```
    CREATE TABLE trajectories (
            seed         INTEGER NOT NULL,
            step         INTEGER NOT NULL,
            reaction_id  INTEGER NOT NULL,
            time         REAL NOT NULL
    );
```

```
    CREATE TABLE factors (
            factor_zero         REAL NOT NULL,
            factor_two          REAL NOT NULL,
            factor_duplicate    REAL NOT NULL)
```

```
    CREATE TABLE initial_state (
            species_id             INTEGER NOT NULL PRIMARY KEY,
            count                  INTEGER NOT NULL
    );

```

## Running NPMC
NPMC is run as follows:

```
NPMC --nano_particle_database=np.sqlite --initial_state_database=initial_state.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200 --dependency_threshold=1
```

- `nano_particle_database`: a sqlite database containing the nano particle data and metadata.
- `initial_state_database` : a sqlite database containing initial state. The simulation trajectories are also written into the database
- `number_of_simulation`: an integer specifying how many simulations to run
- `base_seed`: seeds used are `base_seed, base_seed+1, ..., base_seed+number_of_simulations-1`
- `thread_count`: is how many threads to use.
- `step_cutoff`: how many steps in each simulation

### The Nano particle Database
There are 4 tables in the nano particle database:
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

there are 3 tables in the initial state database:
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

`distance_factor_type` specifies how to compute interaction propensities for two site interactions as a function of distance. Currently the accepted values are `linear` and `inverse_cubic`.

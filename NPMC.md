# NPMC - <span style="color: #0066CC"> Nano Particle Monte Carlo </span>

A three dimensional statistical field theory simulator which supports one and two site interactions. Useful for simulating nanoparticles. Species in this case are dopants to the host matrix. For example, in a nanoparticle composed of a NaYF4 host, any lanthanide such as Yb3+ or Tm3 can be doped onto the Y3+ site. The 4f electrons of these lanthanides give rise to some number of excitation levels. To calculate the rates for the interactions please see [NanoParticleTools](./https://github.com/BlauGroup/NanoParticleTools) as this is a non-trivial process. 

## Sqlite IO

Sqlite is used for input, output, and checkpointing. Before running `NPMC` two necessary .sqlite files must be generated - The Nano Particle Database and State Database. Examples of Python code used to generate these files are available in [Examples](./Examples.html). Below is an outline of each .sqlite file and its necessary tables. **Each .sqlite file must follow this format exactly**. 

### The Nano Particle Database
There are four tables in the nanoparticle database all of which **must be created and filled in by the user**:

- <span style="color:#0066CC"> species </span> : this table consists of one line for each species in the simulation.
    - <span style="color:#006633"> species_id </span>: unique, starts at 0 and must increase in increments of one. Species in this case are dopants to the host matrix.
    - <span style="color:#006633"> degrees_of_freedom </span>: the number of energy levels that lanthanide dopant can access.

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
## The State Database
There are five tables in the initial state database all of which **must be created by the user**: 

- <span style="color:#0066CC"> initial_state </span>: this table represents the initial state of the simulation. **This table must be filled in by the user.**
- <span style="color:#006633"> degree_of_freedom </span>:
the energy level which a dopant atom of the given species is to be initialized to. This is typically 0 for all atoms if initializing a simulation from ground state.
```
CREATE TABLE initial_state (
    site_id            INTEGER NOT NULL PRIMARY KEY,
    degree_of_freedom  INTEGER NOT NULL
);
```

- <span style="color:#0066CC"> trajectories </span>: this table records each interaction executed during the duration of the simulation. For each reaction the seed of the simulation that executed the reaction and corresponding step and time are recorded. 

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
- <span style="color:#0066CC"> factors </span>: **This table must be filled in by the user.**
    - <span style="color:#006633"> distance_factor_type </span>: specifies how to compute interaction propensities for two site interactions as a function of distance. Currently the accepted values are `linear` and `inverse_cubic`.

```
CREATE TABLE factors (
    one_site_interaction_factor      REAL NOT NULL,
    two_site_interaction_factor      REAL NOT NULL,
    interaction_radius_bound         REAL NOT NULL,
    distance_factor_type             TEXT NOT NULL
);
```

- <span style="color:#0066CC"> interrupt_state </span>: during checkpointing, the simulation will fill this table with the final state of the simulation. 

```
CREATE TABLE interrupt_state (
    seed                    INTEGER NOT NULL,
    site_id                 INTEGER NOT NULL,
    degree_of_freedom       INTEGER NOT NULL
        
); 
```

- <span style="color:#0066CC"> interrupt_cutoff </span>: during checkpointing, the simulation will fill in this table.

```
CREATE TABLE interrupt_cutoff (
        seed                    INTEGER NOT NULL,
        step                    INTEGER NOT NULL,
        time                    INTEGER NOT NULL
        
);
```

## Running NPMC
To access the makefile, enter the `NPMC` folder:

```
$ cd `NPMC`
```

Next create an executable with the makefile. The executable will be located in the `NPMC` folder.

```
$ make `NPMC`
```

For further help on the makefile and to view other commands:

```
$ make help
```

`NPMC` requires six input arguments (either `step_cutoff` or `time_cutoff` must be specified): 

- <span style="color:#0066CC"> nano_particle_database </span>: a sqlite database containing the nanoparticle data and metadata.
- <span style="color:#0066CC"> initial_state_database </span> : a sqlite database containing initial state. The simulation trajectories are also written into the database
- <span style="color:#0066CC"> number_of_simulation </span>: an integer specifying how many simulations to run
-  <span style="color:#0066CC"> base_seed </span>: seeds used are `base_seed, base_seed+1, ..., base_seed+number_of_simulations-1`
- <span style="color:#0066CC"> thread_count </span>: how many threads to use.
- <span style="color:#0066CC"> step_cutoff </span>: how many steps in each simulation.
- <span style="color:#0066CC"> time_cutoff </span>: how much time in each simulation [s].

When running `NPMC` ensure that your input file paths are correct considering the executable is inside the `NPMC` folder. Below is an example of how `NPMC` can be run using the input files from [Examples](./Examples.html) (here `step_cutoff` is specified):

```
./NPMC --nano_particle_database=../examples/NPMC/np.sqlite --initial_state_database=../examples/NPMC/initial_state.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200 
```
# GMC 
#### Gillespie Monte Carlo

Implementation of Gillespie's next reaction simulator appropriate for applications in a homogeneous region or where species are well mixed.

## Sqlite IO  

Sqlite is used for input, output, and checkpointing. Before running GMC two necessary .sqlite files must be generated - The Reaction Network Database and State Database. Examples of Python code used to generate these files are available in [Examples](./Examples.html). Below is an outline of each .sqlite file and its necessary tables. **Each .sqlite file must follow this format exactly**. 

### The Reaction Network Database 
There are two tables in the reaction network database both of which **must be created and filled in by the user**:
- <span style="color:#0066CC"> metadata </span> : This table consists of one line for the total number of species and reactions in the simulation.

```
    CREATE TABLE metadata (
            number_of_species   INTEGER NOT NULL,
            number_of_reactions INTEGER NOT NULL
    );
```

- <span style="color:#0066CC"> reactions </span>: This table is how reactions are defined in the simulation. *Only reactions of up to two reactants and products are supported.* Each row in the table represents one reaction with the following attributes. 
    - <span style="color:#006633"> reaction_id </span>: Unique, starts at 0 and must increase in increments of one.
    - <span style="color:#006633"> number_of_reactants\products </span>: Either 0, 1, or 2.
    - <span style="color:#006633"> reactant_1\2 </span>: Unique, positive integer representative of a species. The integer representation of species must begin at 0 and increase in increments of one.
    - <span style="color:#006633"> rate </span>: Rate of the reaction.

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

- <span style="color:#0066CC"> initial_state </span>: This table represents the initial concentration of species. Each row consists of a species_id and corresponding quantity. If there is no row for a species, GMC will initalize its quantity to zero. **This table must be filled in by the user.**

```
    CREATE TABLE initial_state (
            species_id             INTEGER NOT NULL PRIMARY KEY,
            count                  INTEGER NOT NULL
    );
```
- <span style="color:#0066CC"> trajectories </span>: This table records each reaction run during the duration of the simulation. For each reaction the seed of the simulation that executed the reaction and corresponding step and time are recorded. 

```
    CREATE TABLE trajectories (
            seed                INTEGER NOT NULL,
            step                INTEGER NOT NULL,
            time                REAL NOT NULL,
            reaction_id         INTEGER NOT NULL
    );
```
- <span style="color:#0066CC"> factors </span>: This table contains factors that can be used to modify rates of reactions which have zero or two reactants, or have duplicate reactants. **This table must be filled in by the user.**

```
    CREATE TABLE factors (
            factor_zero      REAL NOT NULL,
            factor_two       REAL NOT NULL,
            factor_duplicate REAL NOT NULL
    );
```
- <span style="color:#0066CC"> interrupt_state </span>: During checkpointing, the simulation will fill this table with the final state of the simulation. 

```
    CREATE TABLE interrupt_state (
            seed                    INTEGER NOT NULL,
            species_id              INTEGER NOT NULL,
            count                   INTEGER NOT NULL
            
    );
```
- <span style="color:#0066CC"> interrupt_cutoff </span>: During checkpointing, the simulation will fill in this table.

```
    CREATE TABLE interrupt_cutoff (
            seed                    INTEGER NOT NULL,
            step                    INTEGER NOT NULL,
            time                    INTEGER NOT NULL
            
    );
```
## Running GMC

To run GMC first create an executable with the makefile. 

```
$ make GMC
```
GMC requires six input arguments (either step_cutoff or time_cutoff must be specified): 

- <span style="color:#0066CC"> reaction_database </span>: a sqlite database containing the reaction network and metadata.
- <span style="color:#0066CC"> initial_state_database </span>: a sqlite database containing initial state. The simulation trajectories are also written into the database.
-  <span style="color:#0066CC">number_of_simulation </span>: an integer specifying how many simulations to run.
-  <span style="color:#0066CC">base_seed </span>: seeds used are `base_seed, base_seed+1, ..., base_seed+number_of_simulations-1`.
- <span style="color:#0066CC"> thread_count </span>: how many threads to use.
- <span style="color:#0066CC"> step_cutoff </span>: how many steps in each simulation.
- <span style="color:#0066CC"> time_cutoff </span>: how much time in each simulation [s].

GMC can then be run as follows (here step_cutoff is specified):

```
./GMC --reaction_database=rn.sqlite --initial_state_database=initial_state.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200
```
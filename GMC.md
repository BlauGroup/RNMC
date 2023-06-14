# GMC

Implementation of Gillespie's next reaction simulator appropriate for applications in a homogeneous region or where species are well mixed.

## Sqlite IO

Sqlite is used for input, output, and checkpointing. Before running GMC the necessary .sqlite files must be generated. Examples of Python code used to generate these files are available in [Examples](./Examples.html). Below is an outline of each .sqlite file and its necessary tables. Each .sqlite file **must follow this format exactly**. 

### The Reaction Network Database <span style="color:blue">some *blue* text</span>.
There are 2 tables in the reaction network database:
- `metadata`: This table includes the total number of species and reactions in the simulation.

```
    CREATE TABLE metadata (
            number_of_species   INTEGER NOT NULL,
            number_of_reactions INTEGER NOT NULL
    );
```

- `reactions`: This table is how reactions are defined in the simulation. Only reactions of up to two reactants and products are supported. 
    - `reaction_id`: Unique, starts at 0 and must increase in increments of one
    - `number_of_reactants\products`: Either 0, 1, or 2
    - `reactant_1\2`: Unique, positive integer representative of a species. The integer representation of species must begin at 0 and increase in increments of one
    - `rate`: Rate of the reaction

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
There are 5 tables in the initial state database. 

- `initial_state`: 
- `trajectories`:
- `factors`: The factors can be used to modify rates of reactions which have zero or two reactants, or have duplicate reactants.
- `interrupt_state`: State of the simulation 
- `interrupt_cutoff`: 


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


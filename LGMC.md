# LGMC - <span style="color: #0066CC"> Lattice Gillespie Monte Carlo </span>

Implementation of Gillespie's next reaction simulator which couples a lattice and homogeneous region capable of electrochemical reactions. Either marcus of Butler-Volmer electron transfer theory can be used. The lattice is periodic in the x, y direction and non-periodic in the z direction. The lattice may be static or dynamic.

## Sqlite IO  

Sqlite is used for input, output, and checkpointing. Before running LGMC two necessary .sqlite files must be generated - The Lattice Reaction Network Database and State Database. Examples of Python code used to generate these files are available in [Examples](./Examples.html). Below is an outline of each .sqlite file and its necessary tables. **Each .sqlite file must follow this format exactly**. 

### The Lattice Reaction Network Database 

There are two tables in the lattice reaction network database both of which **must be created and filled in by the user**:
- <span style="color:#0066CC"> metadata </span> : this table consists of one line for the total number of species and reactions in the simulation.

```
     CREATE TABLE metadata (
            number_of_species   INTEGER NOT NULL,
            number_of_reactions INTEGER NOT NULL
    );
```

- <span style="color:#0066CC"> reactions </span>: this table is how reactions are defined in the simulation. *Only reactions of up to two reactants and products are supported.* Each row in the table represents one reaction with the following attributes. 
 - <span style="color:#006633"> reaction_id </span>: unique, starts at 0 and must increase in increments of one.
    - <span style="color:#006633"> number_of_reactants|products </span>: either 0, 1, or 2.
    - <span style="color:#006633"> reactant_1|2 </span>: unique, positive integer representative of a species. The integer representation of species ** must begin at 1 ** and increase in increments of one. ** The integer 0 is reserved to represent an empty site. ** If there is only one reactant|product then set the species to -1.
    - <span style="color:#006633"> phase_reactant|product_1|2 </span>: char representing if the species in the species is in the lattice(spatially resolved), 'L', or solution(homogeneous), 'S', region. If there are not two reactants|products, the phase can be set to 'N'.
    - <span style="color:#006633"> dG </span>: gibbs free energy.
    - <span style="color:#006633"> prefactor </span>: prefactor applied during calculation of reaction rate.
    - <span style="color:#006633"> rate </span>: rate of the reaction.
    - <span style="color:#006633"> electron_tunneling_coefficient </span>: electron tunneling coefficient used for Butler-Volmer and Marcus calculations.
    - <span style="color:#006633"> reorganization_energy </span>: used for Marcus charge transfer.
    - <span style="color:#006633"> charge_transfer_coefficient </span>: charge transfer used for Butler-Volmer and Marcus calculations.
    - <span style="color:#006633"> type </span>: type of reaction:
        - `A`: adsorption
        - `D`: desorption
        - `F`: diffusion (only possible in lattice region)
        - `L`: reaction entirely in the lattice (spatially resolved region)
        - `S`: reaction entirely in the solution (homogeneous region)
        - `O`: oxidation
        - `R`: reduction

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

- <span style="color:#0066CC"> initial_state </span>: this table represents the initial concentration of species. Each row consists of a species_id and corresponding quantity. If there is no row for a species, LGMC will initalize its quantity to zero. **This table must be filled in by the user.**

```
    CREATE TABLE initial_state (
            species_id             INTEGER NOT NULL PRIMARY KEY,
            count                  INTEGER NOT NULL
    );
```
- <span style="color:#0066CC"> trajectories </span>: this table records each reaction run during the duration of the simulation. For each reaction the seed of the simulation that executed the reaction and corresponding step and time are recorded. 
     - <span style="color:#006633"> site_1|2_mapping </span>: single integer representation of the i,j,k values of the lattice site involved in the reaction calculated with the Szudzik algorithm. The code for creating these mappings is shown in [Examples](./Examples.html). The ordering of the sites corresponds to the ordering of the products in the reaction. If one of the products is in the homogeneous region, site_1|2_mapping is equal to -2. If there is only one site involved in the reaction site_1|2_mapping is equal to -3.

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
- <span style="color:#0066CC"> factors </span>: this table contains factors that can be used to modify rates of reactions which have zero or two reactants, or have duplicate reactants. **This table must be filled in by the user.**

```
    CREATE TABLE factors (
            factor_zero      REAL NOT NULL,
            factor_two       REAL NOT NULL,
            factor_duplicate REAL NOT NULL
    );
```
- <span style="color:#0066CC"> interrupt_state </span>: during checkpointing, the simulation will fill this table with the final state of the simulation. This includes both the species in the lattice and homogeneous region.
    - <span style="color:#006633"> site_mapping </span>: szudzik representation of site's i,j,k if on lattice or -2 if in homogeneous region.
    - <span style="color:#006633"> edge </span>: 0 or 1 representative of if the site allows adsorption or desorption reactions.

```
CREATE TABLE interrupt_state (
        seed                    INTEGER NOT NULL,
        species_id              INTEGER NOT NULL,
        quantity                INTEGER NOT NULL,
        site_mapping            INTEGER NOT NULL,
        edge                    INTEGER NOT NULL
        
);
```
- <span style="color:#0066CC"> interrupt_cutoff </span>: During checkpointing, the simulation will fill in this table.

```
    CREATE TABLE interrupt_cutoff (
            seed                    INTEGER NOT NULL,
            step                    INTEGER NOT NULL,
            time                    INTEGER NOT NULL,
            maxk                    INTEGER NOT NULL
            
    );
```
## Running LGMC

To run LGMC first create an executable with the makefile. 

```
$ make LGMC
```
LGMC requires seven input arguments (either step_cutoff or time_cutoff must be specified): 

- <span style="color:#0066CC"> reaction_database </span>: a sqlite database containing the reaction network and metadata.
- <span style="color:#0066CC"> initial_state_database </span>: a sqlite database containing initial state. The simulation trajectories are also written into the database.
-  <span style="color:#0066CC">number_of_simulation </span>: an integer specifying how many simulations to run.
-  <span style="color:#0066CC">base_seed </span>: seeds used are `base_seed, base_seed+1, ..., base_seed+number_of_simulations-1`.
- <span style="color:#0066CC"> thread_count </span>: how many threads to use.
- <span style="color:#0066CC"> step_cutoff </span>: how many steps in each simulation.
- <span style="color:#0066CC"> time_cutoff </span>: how much time in each simulation [s].
- <span style="color:#0066CC"> parameters </span>: name of .txt file containing additional parameters for LGMC
    - <span style="color:#006633"> lattice_constant </span>
    - <span style="color:#006633"> box x upper boundary </span>
    - <span style="color:#006633"> box y upper boundary </span>
    - <span style="color:#006633"> box z upper boundary </span>
    - <span style="color:#006633"> temperature </span>
    - <span style="color:#006633"> Electron free energy </span>
    - <span style="color:#006633"> Is add site (T/F) </span>
    - <span style="color:#006633"> Charge transfer style (M/B) </span>

LGMC can then be run as follows (here step_cutoff is specified):

```
./LGMC --reaction_database=rn.sqlite --initial_state_database=initial_state.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200 --parameters=LGMC_params.txt
```
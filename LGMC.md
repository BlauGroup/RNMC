# LGMC - <span style="color: #0066CC"> Lattice Gillespie Monte Carlo </span>

A kMC implementation coupling a homogeneous (Gillespie-like) region with a lattice, enabling simulations with reactions occurring in multiple phases and capable of electrochemical reactions. Either marcus of Butler-Volmer electron transfer theory (specified at runtime, see 'Running LGMC' below) can be used when calculating rates for electrochemical reactions. The lattice is periodic in the x, y direction and non-periodic in the z direction. The lattice may be static (no sites can be added or deleted) or dynamic (sites are added through adsorption and deleted through desorption reactions).

As a starting point two examples of `LGMC` uses are shown in [Examples](./Examples.html) - `LGMC` with a static lattice to model CO Oxidation on Cu and `LGMC` with a dynamic lattice to simulate the formation and evolution of the solid electrolyte interphase in a lithium-ion battery.

## Sqlite IO  

Sqlite is used for input, output, and checkpointing. Before running LGMC two necessary .sqlite files must be generated - The Lattice Reaction Network Database and State Database. Examples of Python code used to generate these files are available in <a href="{{ site.github.repository_url }}"> examples directory </a>. Below is an outline of each .sqlite file and its necessary tables. **Each .sqlite file must follow this format exactly**. 

### The Lattice Reaction Network Database 

There are two tables in the lattice reaction network database both of which **must be created and filled in by the user**:
- <span style="color:#0066CC"> metadata </span> : this table consists of one line for the total number of species and reactions in the simulation.
<br>
<br>
```
CREATE TABLE metadata (
    number_of_species   INTEGER NOT NULL,
    number_of_reactions INTEGER NOT NULL
);
```
<ul>
<li>
<span style="color:#0066CC"> reactions </span>: this table is how reactions are defined in the simulation. Only reactions of up to two reactants and products are supported. Each row in the table represents one reaction with the following attributes. </li>
    <ul>
    <li> <span style="color:#006633"> reaction_id </span>: unique, starts at 0 and must increase in increments of one. </li>
    <li> <span style="color:#006633"> number_of_reactants/products </span>: either 0, 1, or 2. </li>
    <li> <span style="color:#006633"> reactant_1&#124;2 </span>: unique, positive integer representative of a species. The integer representation of species **must begin at 1** and increase in increments of one. **The integer 0 is reserved to represent an empty site.** If there is only one reactant/product then set the species to -1 </li>
    <li> <span style="color:#006633"> phase_reactant&#124;product_1&#124;2 </span>: char representing if the species in the species is in the lattice(spatially resolved), 'L', or solution(homogeneous), 'S', region. If there are not two reactants/products, the phase can be set to 'N'. </li>
    <li> <span style="color:#006633"> dG </span>: gibbs free energy. </li>
    <li> <span style="color:#006633"> prefactor </span>: prefactor applied during calculation of reaction rate. </li>
    <li> <span style="color:#006633"> rate </span>: rate of the reaction. </li>
    <li> <span style="color:#006633"> electron_tunneling_coefficient </span>: electron tunneling coefficient used for Butler-Volmer and Marcus calculations. </li>
    <li> <span style="color:#006633"> reorganization_energy </span>: used for Marcus charge transfer. </li>
    <li> <span style="color:#006633"> charge_transfer_coefficient </span>: charge transfer used for Butler-Volmer and Marcus calculations. </li>
    <ul>
    <li> <span style="color:#006633"> type </span>: type of reaction: </li>
        <li> `A`: adsorption </li>
        <li> `D`: desorption </li>
        <li> `F`: diffusion (only possible in lattice region) </li>
        <li> `L`: reaction entirely in the lattice (spatially resolved region) </li>
        <li> `S`: reaction entirely in the solution (homogeneous region) </li>
        <li> `O`: oxidation </li>
        <li> `R`: reduction </li>
    </ul> </ul>

<pre><code> CREATE TABLE reactions (
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
</code></pre>
</ul> 

### The State Database 
There are five tables in the initial state database all of which **must be created by the user**: 

- <span style="color:#0066CC"> initial_state </span>: this table represents the initial concentration of species. Each row consists of a species_id and corresponding quantity. If there is no row for a species, LGMC will initalize its quantity to zero. **This table must be filled in by the user.**
<br>
<br>
```
CREATE TABLE initial_state (
        species_id             INTEGER NOT NULL PRIMARY KEY,
        count                  INTEGER NOT NULL
);
```

<ul>
<li> <span style="color:#0066CC"> trajectories </span>: this table records each reaction run during the duration of the simulation. For each reaction the seed of the simulation that executed the reaction and corresponding step and time are recorded. </li>
     <ul>
     <li> <span style="color:#006633"> site_1&#124;2_mapping </span>: single integer representation of the i,j,k values of the lattice site involved in the reaction calculated with the  <a href="./szudzik.md">Szudzik algorithm</a> The ordering of the sites corresponds to the ordering of the products in the reaction. If one of the products is in the homogeneous region, site_1&#124;2_mapping is equal to -2. If there is only one site involved in the reaction site_1&#124;2_mapping is equal to -3. 
     </li>
     </ul>
<pre><code> CREATE TABLE trajectories (
        seed                INTEGER NOT NULL,
        step                INTEGER NOT NULL,
        time                REAL NOT NULL,
        reaction_id         INTEGER NOT NULL,
        site_1_mapping      INTEGER NOT NULL,
        site_2_mapping      INTEGER NOT NULL
);
</code></pre>
</ul>

- <span style="color:#0066CC"> factors </span>: this table contains factors that can be used to modify rates of reactions which have zero or two reactants, or have duplicate reactants. **This table must be filled in by the user.**
<br>
<br>
```
CREATE TABLE factors (
        factor_zero      REAL NOT NULL,
        factor_two       REAL NOT NULL,
        factor_duplicate REAL NOT NULL
);
```

<ul>
<li> <span style="color:#0066CC"> interrupt_state </span>: during checkpointing, the simulation will fill this table with the final state of the simulation. This includes both the species in the lattice and homogeneous region. </li>
    <ul>
    <li> <span style="color:#006633"> site_mapping </span>: szudzik representation of site's i,j,k if on lattice or -2 if in homogeneous region. </li>
    <li> <span style="color:#006633"> edge </span>: 0 or 1 representative of if the site allows adsorption or desorption reactions. </li>
</ul> </ul>

<pre><code> CREATE TABLE interrupt_state (
        seed                    INTEGER NOT NULL,
        species_id              INTEGER NOT NULL,
        quantity                INTEGER NOT NULL,
        site_mapping            INTEGER NOT NULL,
        edge                    INTEGER NOT NULL
        
);
</code></pre>
</ul>

- <span style="color:#0066CC"> interrupt_cutoff </span>: During checkpointing, the simulation will fill in this table.
<br>
<br>
```
CREATE TABLE interrupt_cutoff (
        seed                    INTEGER NOT NULL,
        step                    INTEGER NOT NULL,
        time                    INTEGER NOT NULL,
        maxk                    INTEGER NOT NULL
);
```

## Running LGMC
To access the makefile, enter the LGMC folder:

```
$ cd LGMC
```

Next create an executable with the makefile. The executable will be located in the LGMC folder.

```
$ make LGMC
```

For further help on the makefile and to view other commands:

```
$ make help
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

When running LGMC ensure that your input file paths are correct considering the executable is inside the LGMC folder. Below is an example of how LGMC can be run using the input files from [Examples](./Examples.html) (here `step_cutoff` is specified):

```
./LGMC --lattice_reaction_database=../examples/LGMC/CO_oxidation/rn.sqlite --initial_state_database=../examples/LGMC/CO_oxidation/initial_state.sqlite --number_of_simulations=1000 --base_seed=1000 --thread_count=8 --step_cutoff=200 --parameters=../examples/LGMC/CO_oxidation/LGMC_params.txt
```
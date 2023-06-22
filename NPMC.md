# NPMC - <span style="color: #0066CC"> Nano Particle Monte Carlo </span>

A three dimensional statistical field theory simulator which supports one and two site interactions useful for simulating nanoparticles. Some examples of single site interaction are optical transitions, multiphonon relaxation, or magnetic dipole. Two site interactions represent energy transfer events in which energy is transferred from one species to another.

Species in this case are dopants to the host matrix. For example, in a nanoparticle composed of a NaYF4 host, any lanthanide such as Yb3+ or Tm3 can be doped onto the Y3+ site. The 4f electrons of these lanthanides give rise to some number of excitation levels. Calculating the rates for the interactions is a **non-trivial process**, please refer to [NanoParticleTools](./https://github.com/BlauGroup/NanoParticleTools) for this. 

## Sqlite IO

Sqlite is used for input, output, and checkpointing. Before running `NPMC` two necessary .sqlite files must be generated - The Nano Particle Database and State Database. An example of the Python code used to generate these files is available in <a href="{{ site.github.repository_url }}"> examples directory </a>. Below is an outline of each .sqlite file and its necessary tables. **Each .sqlite file must follow this format exactly**. 

### The Nano Particle Database
There are four tables in the nanoparticle database all of which **must be created and filled in by the user**:

<ul>
<li><span style="color:#0066CC"> species </span> : this table consists of one line for each species in the simulation. </li>
    <ul>
    <li> <span style="color:#006633"> species_id </span>: unique, starts at 0 and must increase in increments of one. Species in this case are dopants to the host matrix. </li>
    <li> <span style="color:#006633"> degrees_of_freedom </span>: the number of energy levels that lanthanide dopant can access. </li>
    </ul>

<pre><code> CREATE TABLE species (
    species_id          INTEGER NOT NULL PRIMARY KEY,
    degrees_of_freedom  INTEGER NOT NULL
);
</code></pre>
</ul>

- <span style="color:#0066CC"> sites </span>: this table initalizes the sites available in the simulation. There are no restrictions in the x, y, z sites table. Atoms may be placed anywhere in space, although it is typically restricted to sites on a lattice. Although that lattice is not fixed, since there are many host materials that are used for upconverting nanoparticles.
```
CREATE TABLE sites (
    site_id             INTEGER NOT NULL PRIMARY KEY,
    x                   REAL NOT NULL,
    y                   REAL NOT NULL,
    z                   REAL NOT NULL,
    species_id          INTEGER NOT NULL
);
```

<ul>
<li> <span style="color:#0066CC"> interactions </span>: Interactions are the energy transitions which take species from one energy level to another. </li> 
    <ul>
    <li> <span style="color:#006633"> interaction_id </span>: unique, index which monotonically increases starting from 0. </li>
    <li> <span style="color:#006633"> number_of_sites </span>: number of sites which participate in the event, either 1 or 2. </li>
    <li> <span style="color:#006633"> species_id_1\|2 </span>: species_id corresponding to definitions provided in species table. If a single site interaction, species_id_2 should be -1. </li>
    <li> <span style="color:#006633"> left_state_1\|2 </span>: corresponds to the initial energy level of a species (analogous to left side of a reaction). If a single site interaction, left_state_2 should be -1. </li>
    <li> <span style="color:#006633"> right_state_1\|2 </span>: corresponds to the final energy level of a species (analogous to right side of a reaction). If a single site interaction, right_state_2 should be -1. </li>
    <li> <span style="color:#006633"> rate </span>: rate for the energy transition event, please refer to [NanoParticleTools](./https://github.com/BlauGroup/NanoParticleTools). </li> </ul>

<pre><code> CREATE TABLE interactions (
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
</code></pre>
</ul>

- <span style="color:#0066CC"> metadata </span>: this table consists of one line for the total number of species, sites, and interactions in the simulation.
```
CREATE TABLE metadata (
    number_of_species                   INTEGER NOT NULL,
    number_of_sites                     INTEGER NOT NULL,
    number_of_interactions              INTEGER NOT NULL
);
```
## The State Database
There are five tables in the initial state database all of which **must be created by the user**: 

<ul>
<li> <span style="color:#0066CC"> initial_state </span>: this table represents the initial state of the simulation. <b>This table must be filled in by the user.</b> </li>
   <ul>
    <li> <span style="color:#006633"> degree_of_freedom </span>: </li>
    the energy level which a dopant atom of the given species is to be initialized to. This is typically 0 for all atoms if initializing a simulation from ground state. </ul>
<pre><code> CREATE TABLE initial_state (
    site_id            INTEGER NOT NULL PRIMARY KEY,
    degree_of_freedom  INTEGER NOT NULL
);
</code></pre>
</ul>

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

<ul>
<li> <span style="color:#0066CC"> factors </span>: this table contains factors that can be used to modify the rates of interactions. <b> This table must be filled in by the user.</b> </li>
    <ul> <li> <span style="color:#006633"> distance_factor_type </span>: specifies how to compute interaction propensities for two site interactions as a function of distance. Currently the accepted values are `linear` and `inverse_cubic`. </li>

<pre><code> CREATE TABLE factors (
    one_site_interaction_factor      REAL NOT NULL,
    two_site_interaction_factor      REAL NOT NULL,
    interaction_radius_bound         REAL NOT NULL,
    distance_factor_type             TEXT NOT NULL
);
</code></pre>
</ul>

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
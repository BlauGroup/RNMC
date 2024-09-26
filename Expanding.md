# Adding Simulators to RNMC

`RNMC` is designed to be modular. Individual types of simulations all draw from a small `core` library, but most of their logic is independent. This means that, to develop a new type of simulation, you do not need to have a deep understanding of every component of `RNMC`. As long as you understand `core` and the general structure of a new simulation type, you're free to do what you want!

Here, we describe the structure of `RNMC`, explaining the components of the `core` library and what goes into a simulator in order to help facilitate the development of new simulators.

## The Core Library

The `core` library contains resources that are shared by just about every simulator in `RNMC`. New simulators may require small modifications to `core`, but it also may be possible for you to develop a new simulator without changing anything!

`core` contains:

- `RNMC_types.h`: Commonly used datatypes for simulations.
- `sql.h`: An interface to read from and write to SQLite3 database files.
- `sql_types.h` / `sql_types.cpp`: Specific datatypes for interacting with types of tables and data that are commonly used by `RNMC` simulators. Examples include the `initial_state` table (defining what the initial simulation looks like, in terms of numbers of each species) and `metadata` table (defining how many different species and reactions are in the simulation).
- `sampler.h`: A wrapper over `GSL` random number generators, used for random sampling.
- `queues.h`: Queues for parallel operation, ensuring that different processes do not attempt to run simulations on the same seed and that history data is safely added to the trajectory database.
- `simulation.h` / `simulation.cpp`: Defines a core interface to execute steps with either a time or number-of-steps cutoff. Simulations on a lattice behave somewhat differently and are defined in `lattice_simulation.h` and `lattice_simultion.cpp`.
- `simulator_payload.h`: A layer above the simulation interface that initializes simulations, executes them, and then logs their outputs to a history queue to be inserted into a database.
- `dispatcher.h` / `dispatcher.cpp`: Defines a `Dispatcher` that manages different simulating processes.

The `core` folder also defines the high-level interfaces for each type of simulator. For instance, `nano_particle_simulation.h`/`nano_particle_simulation.cpp` defines the interface for `NPMC`.

## Components of a Simulator

Different simulators may have different levels of complexity and require different components. Generally, though, `RNMC` simulations have similar components. These include:
- Additional types for SQLite3 input and output
- Additional non-SQL datatypes
- "Solvers", which use the Sampler interface defined in `core/sampler.h` to select events (see *e.g.*, `GMC/linear_solver.h` and `GMC/linear_solver.cpp` for a simple example).
- A "reaction network" class. These classes are expected to manage:
    - Initializing the simulation state by reading from an appropriate SQLite table
    - Updating the state and event propensities based on an event selected by a Solver through the simulation interface
    - Checkpointing simulations
    - Creating history elements to write to a SQLite database
- A command-line interface for simulations. In general, we prefer to define simulation parameters using command-line arguments, but, depending on the complexity of the simulation inputs, reading parameters from a configuration file may also be appropriate.
- A makefile for compiling the simulation executable.

Other than the main simulation interface mentioned above (which is housed in `core`) and testing code (housed in `tests`), all application-specific code should be placed in a separate directory (*e.g.*, `NPMC` or `GMC`). For a more detailed understanding of how `RNMC` is organized, we encourage you to select an existing type of simulation (*e.g.*, `GMC`) and read through the source code.
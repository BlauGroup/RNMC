# Introduction 
Reaction Network Monte Carlo (RNMC) is a collection of programs for kinetic Monte Carlo (kMC) simulation of physical systems heavily based on SPPARKS. RNMC is designed to run large numbers of simulations of a fixed system in parallel. The project consists of three kMC modules for different domains which use `core` code for shared processes, for example IO, threading logic and model independent simulation logic.

## Three Modules
- `GMC` - [Gillespie Monte Carlo](./GMC.md):
Implementation of Gillespie's next reaction simulator. GMC is able to run simulations of reaction networks with hundreds of millions of reactions.
- `NPMC` - [Nano Particle Monte Carlo](./NPMC.md): 
A three dimensional statistical field theory simulator which supports one- and two-site interactions. Useful for simulating nano particles.
- `LGMC` - [Lattice Gillespie Monte Carlo](./NPMC.md): A kMC implementation coupling a homogeneous (Gillespie-like) region with a lattice, enabling simulations with reactions occurring in multiple phases and capable of electrochemical reactions.

# Introduction 
Reaction Network Monte Carlo (RNMC) is a collection of programs for Monte Carlo simulation of statistical mechanical systems heavily inspired by SPPARKS. RNMC is designed to run large numbers of simulations of a fixed system in parallel. The project consists of three monte carlo algorithms for different domains which use `core` code for shared processes, for example IO, threading logic and model independent simulation logic.

## Three Modules
- `GMC` [Gillespie Monte Carlo](./GMC.md): Implementation of Gillespie's next reaction simulator. GMC is able to run simulations of reaction networks with hundreds of millions of reactions, even when the number of species is small.
- `NPMC` [Nano Particle Monte Carlo](./NPMC.md): A three dimensional statistical field theory simulator which supports one and two site interactions. Useful for simulating nano particles.
- `LGMC` [Lattice Gillespie Monte Carlo](./NPMC.md):  Implementation of Gillespie's next reaction simulator which couples a lattice and homogeneous region capable of electrochemical reactions.

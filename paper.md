--
title: 'RNMC: kinetic Monte Carlo implementations for complex reaction networks'
tags:
  - Python
  - C++
  - chemical dynamics
  - kinetic Monte Carlo
authors:
  - name: Laura Zichi
    equal-contrib: true
    corresponding: false
    affiliation: 1,2
  - name: Daniel Barter
    corresponding: false
    equal-contrib: true
    affiliation: 3
  - name: Eric Sivonxay
    corresponding: false
    equal-contrib: true
    affiliation: 3
  - name: Evan Spotte-Smith
    corresponding: false  
    affiliation: 1,4
  - name: Rohith Srinivaas M
    corresponding: false  
    affiliation: 1,4
  - name: Kristin Persson
    corresponding: true  
    affiliation: 4, 5
  - name: Sam Blau
    corresponding: true  
    affiliation: 3

affiliations:
 - name: Materials Science Division, Lawrence Berkeley National Laboratory, Berkeley, CA, USA 94720
   index: 1
 - name: Department of Physics, University of Michigan - Ann Arbor, Ann Arbor, MI, USA 48109
   index: 2
 - name: Energy Storage and Distributed Resources, Lawrence Berkeley National Laboratory, Berkeley, CA USA 94720
   index: 3
 - name: Department of Materials Science and Engineering, University of California - Berkeley, CA, USA 94720
   index: 4
 - name: Molecular Foundry, Lawrence Berkeley National Laboratory, Berkeley, CA, USA 94720
   index: 5
 
date: 30 June 2023
bibliography: paper.bib


# Summary

Macroscopic chemical and physical phenomena are driven by microscopic interactions at the atomic and molecular scales.
In order to capture complex processes with high fidelity, simulation methods that bridge disparate time and length scales are needed.
While techniques like molecular dynamics and ab initio simulations capture dynamics and reactivity at high resolution, they cannot be used beyond relatively short time scales (picoseconds to nanoseconds).
The kinetic Monte Carlo (kMC) algorithm bridges length and time scales across several orders of magnitude, while retaining relevant microscopic resoluton, rendering it a powerful and flexible tool.

`RNMC` is an easy-to-use, high-performance kMC simulation framework that enables modeling of complex systems that are inaccessible with current methods.
It consists of a core module defining the common features of kMC algorithms including an implementation of the Gillespie algorithm [@gillespie1977exact], input/output operations based on SQLite databases, threading logic for parallel execution, and dependency graphs for efficient propensity updates.
In addition, there are currently three modules defining kMC implementations for different domains.
The `GMC` (Gillespie Monte Carlo) module enables simulations of reaction networks in a homogeneous (well-mixed) environment.
The `NPMC` (Nano Particle Monte Carlo) module enables simulation of dynamics in nanoparticles with 3D statistical field theory and supports one- and two-site interactions.
Finally, the `LGMC` (Lattice Gillespie Monte Carlo) module is designed for simulations of multi-phase systems (especially at solid-fluid interfaces) where chemical and electrochemical reactions can occur between a lattice region and a homogeneous region.
We have designed `RNMC` to be easily extensible, enabling users to add additional kMC modules for other diverse chemical and physical systems.


# Statement of need

`RNMC` was designed to address several critical gaps in available open-source kMC simulations (e.g. the Stochastic Parallel PARticle Kinetic Simulator or `SPPARKS`[@garcia2009crossing] and `kmos`).[@hoffmann2014kmos]
In particular, previously developed kMC codes are not well designed for large reaction networks and typically aim to simulate a small set of relatively simple systems (e.g. well-mixed fluids, crystal lattices, or static two-dimensional surfaces).
The module `GMC` can easily scale to hundreds of millions of reactions, as we have demonstrated in prior works on Li-ion and Mg-ion batteries.[@spotte2022toward; @barter2023predictive; @spotte2023chemical].
`NPMC` provides a 3D kMC simulation on a lattice that is specifically designed to model complex interactions in nanocrystals.
Finally, most available kMC implementations cannot easily simulate multi-phase systems or electrochemical processes.
The `LGMC` module overcomes these restrictions and can simulate both coupled reactions between homogeneous and lattice regions with electrochemical reactions based on Marcus[@marcus1965theory] or Butler-Volmer kinetics.[@newman2021electrochemical]
`RNMC` offers a high-performance kMC implementation and specifically enables exciting opportunities for dynamical modeling of previously inaccessible domains.


# Acknowledgements

This project was intellectually led by the Laboratory Directed Research and Development Program of Lawrence Berkeley National Laboratory under U.S. Department of Energy Contract No. DE-AC02-05CH11231.
E.W.C.S.-S. was supported by the Kavli Energy NanoScience Institute Philomathia Graduate Student Fellowship.
Additional support came from the Joint Center for Energy Storage Research (JCESR), an Energy Innovation Hub funded by the U.S. Department of Energy, Office of Science, Basic Energy Sciences.
This code was developed and tested using computational resources provided by the National Energy Research Scientific Computing Center (NERSC), a U.S. Department of Energy Office of Science User Facility under Contract No. DE-AC02-05CH11231, the Eagle and Swift HPC systems at the National Renewable Energy Laboratory (NREL), and the Lawrencium HPC cluster at Lawrence Berkeley National Laboratory.

# References

---
title: >-
    RNMC: kinetic Monte Carlo implementations for complex reaction networks
tags:
  - C++
  - chemical dynamics
  - kinetic Monte Carlo
  - nanoparticle
  - electrochemistry
  - Gillespie
authors:
  - name: Laura Zichi
    orcid: 0000-0003-3897-3097
    equal-contrib: true
    corresponding: false
    affiliation: "1, 2"
  - name: Daniel Barter
    orcid: 0000-0002-6408-1255
    corresponding: false
    equal-contrib: true
    affiliation: 3
  - name: Eric Sivonxay
    orcid: 0000-0002-6408-1255
    corresponding: false
    equal-contrib: true
    affiliation: 3
  - name: Evan Walter Clark Spotte-Smith
    orcid: 0000-0003-1554-197X
    corresponding: false  
    affiliation: "1, 4"
  - name: Rohith Srinivaas Mohanakrishnan
    corresponding: false  
    affiliation: "1, 4"
  - name: Emory M. Chan
    corresponding: false
    affiliation: 5
  - name: Kristin Aslaug Persson
    corresponding: true  
    affiliation: "4, 5"
  - name: Samuel M. Blau
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
date: 14 August 2024
bibliography: paper.bib
---

# Summary
Macroscopic chemical and physical phenomena are driven by microscopic interactions at the atomic and molecular scales.
In order to capture complex processes with high fidelity, simulation methods that bridge disparate time and length scales are needed.
While techniques like molecular dynamics and *ab initio* simulations capture dynamics and reactivity at high resolution, they cannot be used beyond relatively small length (hundreds to thousands of atoms) and time scales (picoseconds to microseconds).
Kinetic Monte Carlo (kMC) approaches overcome these limitations to bridge length and time scales across several orders of magnitude while retaining relevant microscopic resolution, making it a powerful and flexible tool.


Here, we present `RNMC`, an easy-to-use, modular, high-performance kMC simulation framework that enables modeling of complex systems.
`RNMC` consists of a core module defining the common features of kMC algorithms, including an implementation of the Gillespie algorithm [@gillespie1977exact], input/output operations leveraging SQLite databases, threading logic for parallel execution, and dependency graphs for efficient event propensity updates.
In addition, there are currently three modules defining kMC implementations for different types of applications.
The `GMC` (Gillespie Monte Carlo) module enables simulations of reaction networks in a homogeneous (well-mixed) environment.
`GMC` is a basic tool that is appropriate for general simulations of solution-phase chemistry.
The `NPMC` (NanoParticle Monte Carlo) module enables simulation of dynamics in nanoparticles with 3D statistical field theory and supports one- and two-site interactions.
Finally, the `LGMC` (Lattice Gillespie Monte Carlo) module is designed for simulations of multi-phase systems (especially at solid-fluid interfaces) where chemical and electrochemical reactions can occur between a lattice region and a homogeneous region.
We have designed `RNMC` to be easily extensible, enabling users to add additional kMC modules for other diverse chemical and physical systems.

# Statement of need

Three are many existing kMC implementations, including several open source examples (e.g. the Stochastic Parallel PARticle Kinetic Simulator or `SPPARKS` [@garcia2009crossing] and `kmos` [@hoffmann2014kmos]).
`RNMC` began as a fork of SPPARKS but differs in several important ways.
First, because `RNMC` uses the widely supported SQLite database engine for simulation inputs and outputs, it facilitates the automation of simulations.
Second, `RNMC` has a focus on modularity; it is designed such that users can quickly develop new types of kMC simulations using a common core library.
 
The simulation modules already implemented in `RNMC` provide unique capabilities that are not widely available in other open source codes.
`NPMC` is specifically designed for 3D simulations of the complex photophysical interaction networks in nanocrystals [@teitelboim2019energy], particularly multi-domain heterostructures whose optical properties cannot be calculated deterministically [@skripka2023NL].
`NPMC` can be used to simulate energy transfer interactions between dopants in nanoparticles, their radiative transitions, and nonlinear processes such as upconversion [@chan2015combinatorial] and photon avalanching [@skripka2023NL].  
`LGMC` is also somewhat unique in that it can simulate multi-phase systems and electrochemical processes.
Simulations using `LGMC` can include a lattice region and a homogeneous solution region which can interact *via* interfacial reactions.
Electrochemcial reactions can be treated using Marcus theory [@marcus1965theory] or Butler-Volmer kinetics [@newman2021electrochemical].
Because it allows for a dynamic lattice region, `LGMC` is also appropriate for simulations of nucleation and growth, dissolution, precipitation, and related phenomena.

We have already used the `GMC` module in a number of prior works in applications related to Li-ion and Mg-ion batteries [@spotte2022toward; @barter2023predictive; @spotte2023chemical]. We note that these simulations included tens of millions of reactions, demonstrating that `RNMC` is able to scale to large and complex reaction networks. In addition, we have used `NPMC` to perform Bayesian optimization of upconverting nanoparticles [@xia2023accelerating].

# Acknowledgements

This project was intellectually led by the Laboratory Directed Research and Development Program of Lawrence Berkeley National Laboratory under U.S. Department of Energy Contract No. DE-AC02-05CH11231.
L.Z. was supported in part by the U.S. Department of Energy, Office of Science, Office of Workforce Development for Teachers and Scientists (WDTS) under the Science Undergraduate Laboratory Internships Program (SULI).
E.W.C.S.-S. was supported by the Kavli Energy NanoScience Institute Philomathia Graduate Student Fellowship.
Work at the Molecular Foundry (E.M.C., K.A.P) was supported by the Office of Science, Office of Basic Energy Sciences, of the U.S. Department of Energy under Contract No. DE-AC02-05CH11231.
Additional support came from the Joint Center for Energy Storage Research (JCESR), an Energy Innovation Hub funded by the U.S. Department of Energy, Office of Science, Basic Energy Sciences.
This code was developed and tested using computational resources provided by the National Energy Research Scientific Computing Center (NERSC), a U.S. Department of Energy Office of Science User Facility under Contract No. DE-AC02-05CH11231, the Eagle and Swift HPC systems at the National Renewable Energy Laboratory (NREL), and the Lawrencium HPC cluster at Lawrence Berkeley National Laboratory.

# References


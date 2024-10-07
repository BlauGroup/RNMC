<img src="./logo.png">

Reaction Network Monte Carlo (`RNMC`) is a collection of programs for Monte Carlo simulation of statistical mechanical systems heavily inspired by [SPPARKS](https://spparks.sandia.gov/). `RNMC` is designed to run large numbers of simulations of a fixed system in parallel. The project currently consists of four parts:
- `core` : Core code shared by all simulators, for example IO, threading logic, and model independent simulation logic.
- `GMC` : Implementation of Gillespie's next reaction simulator. GMC is able to run simulations of reaction networks with hundreds of millions of reactions, even when the number of species is small.
- `NPMC` : A 3D statistical field theory simulator which supports one and two site interactions. Useful for simulating nano particles.
- `LGMC` : A simulator that can include a static or dynamic lattice region and a homogeneous region and allows electrochemical as well as chemical reactions. Suitable for multi-phase simulations (e.g., heterogeneous catalysis, electrochemical plating or stripping).

Examples of research projects using `RNMC`:
- Spotte-Smith, Kam, et al., *ACS Energy Lett.* **7**(4), 1446–1453 (2022). [DOI: 10.1021/acsenergylett.2c00517](https://doi.org/10.1021/acsenergylett.2c00517)
- Barter, Spotte-Smith, et al., *Dig. Disc.* **2**, 123-137 (2023). [DOI: 10.1039/D2DD00117A](https://doi.org/10.1039/D2DD00117A)
- Spotte-Smith et al., *J. Am. Chem. Soc.* **145**(2), 12181–12192 (2023). [DOI: 10.1021/jacs.3c02222](https://doi.org/10.1021/jacs.3c02222)
- Xia et al., *Nano Lett.* **23**(23), 11129–11136 (2023). [DOI: 10.1021/acs.nanolett.3c03568](https://doi.org/10.1021/acs.nanolett.3c03568)

## Documentation

The complete documentation for `RNMC` can be found [here](https://blaugroup.github.io/RNMC/). This includes a guide to installation, setting up simulations, testing, and more.

## Dependencies

RNMC depends on [GSL](https://www.gnu.org/software/gsl/) for pseudo random number generation and [sqlite](https://www.sqlite.org/index.html) for the database interfaces.

## Contributing

We welcome community contributions to RNMC! Please submit code contributions as [pull requests (PRs)](https://github.com/BlauGroup/RNMC/pulls). Work-in-progress PRs are encouraged; just include "\[WIP\]" in the title of your PR.

If you don't know where to start, you can check out the open [issues](https://github.com/BlauGroup/RNMC/issues) or reach out to [Sam Blau](mailto:smblau@lbl.gov) and describe what ideas you may have in mind.

Have a problem installing or running RNMC, or have an idea for how the code could be better? Need support in using RNMC? Please open a new issue, and label the issue (with *e.g.*, “bug”, “enhancement”, or “question”) so that we can triage appropriately. Issues are preferred over e-mails or other private communications because multiple users might encounter the same problem.

A more complete guide to contributing can be found in our [online documentation](https://blaugroup.github.io/RNMC/Contributors.html).
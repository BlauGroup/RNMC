# Examples

## GMC

## NPMC

## LGMC
Two examples of applications of LGMC, CO Oxidation on Cu with a static lattice and a simple reaction network for the formation and evolution of the solid electrolyte interphase (SEI) in a lithium-ion battery, are provided. The Pyhton code for generating the .sqlite files is in the `analysis` folder of the github. 

### [CO Oxidation on Cu](./CO_oxidation.html)
This example simulates electrocatalytic CO oxidation on Cu. We employ a static lattice with 50 sites in each of the x and y dimensions and 2 sites in the z dimension to represent a Cu surface as the catalysis. 

All initial sites are empty. The initial state of the solution in contact with this surface consists of 2,500 CO molecules and 15,000 H_2O molecules. The reduction and oxidation rates use an electron free energy of -0.5 eV. The system is modeled at 300K. 

The allowed reactions and associated rates are as follows:

| Reaction                                 | Prefactor | Rate $s^{-1}$ |
|------------------------------------------|-----------|---------------|
| $CO_{soln}$ + * &rarr; $CO^*$            | 1         | $10^5$        |
| * + $H_2O_{(soln)}$ &rarr; $H_2O^*$      | 1         | $10^2$        |
| $H_2O^*$ &rarr; * + $H_2O_{(soln)}$      | 1         | $10^2$        |
| $H_2O^*$   &rarr;  $OH^*$ + $e^-$        | 0.02      | 0             |
| $OH^*$ + $e^-$ &rarr; $H_2O^*$           | $10^4$    | 0             |
| $CO^*$ + $OH^*$ &rarr; $CO_2^*$ + $e^-$  | 0.8432    | 0             |
| $CO_2^*$ &rarr; * + $CO_2_{(soln)}$      | 1         | $10^4$        |
| $CO^*$ + * &rarr; * + $CO^*$             | 1         | 1             |
| $OH^*$ + * &rarr; * + $OH^*$             | 1         | 1             |

From this initial state, 200,000 steps of our kMC are run. As the simulation proceeds, CO and H~2~O rapidly adsorbs onto the lattice. The H~2~O on the lattice then oxidizes to form OH which reacts with the CO on the lattice to form CO~2~. This CO ~2~ lattice product then desorbs into the solution. The results of the simulation are shown as the occupancy of the lattice sites.

<figure>
    <img src="valid.png"
         alt="Results of CO Oxidation on Cu">
    <figcaption>Surface occupancy for CO and OH where empty represents a site that contains neither CO or OH. </figcaption>
</figure>
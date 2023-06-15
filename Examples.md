# Examples

## GMC

## NPMC

## LGMC
Two examples of applications of LGMC, [CO Oxidation on Cu](./CO_oxidation.html) with a static lattice and a simple reaction network for the formation and evolution of the [solid electrolyte interphase (SEI)](./SEI.md) in a lithium-ion battery, are provided. The Pyhton code for generating the .sqlite files is in the `analysis` folder of the github. 

### [CO Oxidation on Cu](./CO_oxidation.html)
This example simulates electrocatalytic CO oxidation on Cu. We employ a static lattice with 50 sites in each of the x and y dimensions and 2 sites in the z dimension to represent a Cu surface as the catalysis. 

All initial sites are empty. The initial state of the solution in contact with this surface consists of 2,500 CO molecules and 15,000 H_2O molecules. The reduction and oxidation rates use an electron free energy of -0.5 eV. The system is modeled at 300K. 

The allowed reactions and associated rates are as follows:

| Reaction                                                                | Prefactor | Rate ($s^{-1}$) |
|-------------------------------------------------------------------------|-----------|-----------------|
| $\text{CO}_{(soln)} + * \xrightarrow[]{} \text{CO}^*$                   | 1         | $10^5$          |
| $* + \text{H}_2\text{O}_{(soln)} \xrightarrow[]{} \text{H}_2\text{O}^*$ | 1         | $10^2$          |
| $\text{H}_2\text{O}^* \xrightarrow[]{} * + \text{H}_2\text{O}_{(soln)}$ | 1         | $10^2$          |
| $\text{H}_2\text{O}^* \xrightarrow[]{} \text{OH}^* + e^-$               | 0.02      | 0               |
| $\text{OH}^* + e^- \xrightarrow[]{} \text{H}_2\text{O}^*$               | $10^4$    | 0               |
| $\text{CO}^* + \text{OH}^* \xrightarrow[]{} \text{CO}_2^* + e^-$        | 0.8432    | 0               |
| $\text{CO}_2^* \xrightarrow[]{} * + {\text{CO}_2}_{(soln)}$             | 1         | $10^4$          |
| $\text{CO}^* + * \xrightarrow[]{} * + \text{CO}^*$                      | 1         | 1               |
| $\text{OH}^* + * \xrightarrow[]{} * + \text{OH}^*$                      | 1         | 1               |

From this initial state, 200,000 steps of our kMC are run. As the simulation proceeds, CO and H~2~O rapidly adsorbs onto the lattice. The H~2~O on the lattice then oxidizes to form OH which reacts with the CO on the lattice to form CO~2. This CO~2 lattice product then desorbs into the solution. The results of the simulation are shown as the occupancy of the lattice sites.

![Results of CO Oxidation on Cu](valid.ping)

<figure>
    <img src="valid.jpg"
         alt="Results of CO Oxidation on Cu">
    <figcaption>Surface occupancy for CO and OH where empty represents a site that contains neither CO or OH. </figcaption>
</figure>
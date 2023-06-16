# Examples

## GMC

## NPMC

## LGMC
Two examples of applications of LGMC, CO Oxidation on Cu with a static lattice and a simple reaction network for the formation and evolution of the solid electrolyte interphase (SEI) in a lithium-ion battery, are provided. The Python code for generating the .sqlite files is in the `examples` folder of the github. 

### CO Oxidation on Cu
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

## SEI
This example simulates a simplified formation and evolution of the solid electrolyte interphase in a lithium-ion battery.


### Species
EC0, EC-, LiEC+, LiEC0, LiEC_RO0, LiEC_RO-,
LiCO3-, Li2CO30, LEDC0, LEDC-, LEDC_minus_Li-, LEDC_plus_Li+, LEDC_plus_Li0, C2H40

### Overview:
The molecular thermodynamics are in "test_species_thermo.json". Non-electrochemical reactions are in "test_energy_barriers.json"; inner reorganization energies for reduction reactions are in "test_lambda_inner.json". For outer reorganization energies, we can just use a constant value for now 0.32 eV should be reasonable. For the hopping reactions (in "test_hopping_reactions.json"), I defined the barrier based on what Li was hopping from. For LiEC+, for instance, the barrier is 0.27 eV (based on the residence time of Li+ with EC solvation shells), but for LEDC0 it's 0.64 eV (based on values taken from the literature).

There are a couple of considerations not included in these files:

1. Phase: I didn't mark what phase each reaction can occur in. I think it's probably fine to allow all reactions to occur in either phase, though certain reactions will be unlikely to happen in one phase or another because of adsorption/desorption.
2. Adsorption/desorption: I didn't include adsorption or desorption reactions. In terms of species that we want to be able to adsorp, the most important are the solid products (LiCO3-, Li2CO30, LEDC0, LEDC-, LEDC_minus_Li-, LEDC_plus_Li+, LEDC_plus_Li0). These should adsorb very favorably (adsorption should be fast; desorption should be slow). Technically any of the EC-like species (EC0, EC-, LiEC+, LIEC0, LiEC_RO0, LiEC_RO-) should also be able to adsorb, but to keep the model simple we could just not allow that. In terms of desorption, the only thing that should desorb fast is C2H4. I don't have a good idea for rates for adsorption/desorption reactions; you might need to play around with it.
3. Lattice diffusion: Through the hopping reactions, we can capture lattice diffusion of Li+ (which is the most important). I didn't look up diffusion of CO3^2- or EDC^2- anions. To keep this model as simple as possible, we could just not include any lattice diffusion. Or if you can find some reasonable barrier/rate from the literature, that works too. I'm easy here.
4. Sizes: 
    1: C2H40, LiCO3-, Li2CO30
    2: EC0, EC-, LiEC+, LiEC0, LiEC_RO0, LiEC_RO-
    3: LEDC0, LEDC-, LEDC_minus_Li-, LEDC_plus_Li+, LEDC_plus_Li0
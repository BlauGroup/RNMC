/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#include "nano_solver.h"

NanoSolver::NanoSolver(
    unsigned long int seed,
    std::vector<NanoReaction> &&current_reactions) : sampler(Sampler(seed)),
                                                     cumulative_propensities(current_reactions.size()),
                                                     number_of_active_indices(0),
                                                     propensity_sum(0.0),
                                                     // if this move isn't here, the semantics is that initial
                                                     // propensities gets moved into a stack variable for the function
                                                     // call and that stack variable is copied into the object.
                                                     current_reactions(std::move(current_reactions))
{
    update();

} // NanoSolver()

/* ---------------------------------------------------------------------- */

NanoSolver::NanoSolver(
    unsigned long int seed,
    std::vector<NanoReaction> &current_reactions) : sampler(Sampler(seed)),
                                                    cumulative_propensities(current_reactions.size()),
                                                    number_of_active_indices(0),
                                                    propensity_sum(0.0),
                                                    current_reactions(current_reactions)
{
    update();

} // NanoSolver()

/* ---------------------------------------------------------------------- */

void NanoSolver::update()
{
    cumulative_propensities.resize(current_reactions.size());
    number_of_active_indices = cumulative_propensities.size();
    // std::cerr << "Propensities are: ";
    if (number_of_active_indices > 0)
    {
        cumulative_propensities[0] = current_reactions[0].rate;
        // std::cerr << current_reactions[0].rate << ", ";
    }

    for (unsigned int i = 1; i < current_reactions.size(); i++)
    {
        cumulative_propensities[i] = cumulative_propensities[i - 1] + current_reactions[i].rate;
        //   std::cerr << current_reactions[i].rate << ", ";
    }
    // std::cerr << "\n";

    propensity_sum = cumulative_propensities[cumulative_propensities.size() - 1];

} // update()

/* ---------------------------------------------------------------------- */

std::optional<Event> NanoSolver::event()
{
    if (number_of_active_indices == 0)
    {
        propensity_sum = 0.0;
        return std::optional<Event>();
    }

    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;

    unsigned long m;

    // This is a linear search. In the future, use a binary search
    for (m = 0; m < current_reactions.size(); m++)
    {
        if (cumulative_propensities[m] > fraction)
            break;
    }
    // std::cerr << "Propensity sum: " << propensity_sum << "\n";

    double dt = -std::log(r2) / propensity_sum;
    if (m < current_reactions.size())
        return std::optional<Event>(Event{.index = m, .dt = dt});
    else
        return std::optional<Event>(Event{.index = last_non_zero_event, .dt = dt});

} // event()

/* ---------------------------------------------------------------------- */

double NanoSolver::get_propensity(int index)
{
    return current_reactions[index].rate;

} // get_propensity()

/* ---------------------------------------------------------------------- */

double NanoSolver::get_propensity_sum()
{
    return propensity_sum;

} // get_propensity_sum()
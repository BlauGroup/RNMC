/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#include <vector>
#include <optional>
#include <cmath>
#include <map>

#include "linear_solver.h"

// LinearSolver can opperate directly on the passed propensities using a move
LinearSolver::LinearSolver(
    unsigned long int seed,
    std::vector<double> &&initial_propensities) : sampler(Sampler(seed)),
                                                  // if this move isn't here, the semantics is that initial
                                                  // propensities gets moved into a stack variable for the function
                                                  // call and that stack variable is copied into the object.
                                                  propensities(std::move(initial_propensities)),
                                                  number_of_active_indices(0),
                                                  propensity_sum(0.0)
{
    for (unsigned long i = 0; i < propensities.size(); i++)
    {
        propensity_sum += propensities[i];
        if (propensities[i] > 0)
        {
            number_of_active_indices += 1;
            last_non_zero_event = i;
        }
    }
} // LinearSolver()

/*---------------------------------------------------------------------------*/

LinearSolver::LinearSolver(
    unsigned long int seed,
    std::vector<double> &initial_propensities) : sampler(Sampler(seed)),
                                                 propensities(initial_propensities),
                                                 number_of_active_indices(0),
                                                 propensity_sum(0.0)
{

    for (unsigned long i = 0; i < propensities.size(); i++)
    {
        propensity_sum += propensities[i];
        if (propensities[i] > 0)
        {
            number_of_active_indices += 1;
            last_non_zero_event = i;
        }
    }
} // LinearSolver()

/*---------------------------------------------------------------------------*/

void LinearSolver::update(Update update)
{

    if (propensities[update.index] > 0.0)
        number_of_active_indices--;

    if (update.propensity > 0.0)
    {
        number_of_active_indices++;
        if (update.index > last_non_zero_event)
            last_non_zero_event = update.index;
    }

    propensity_sum -= propensities[update.index];
    propensity_sum += update.propensity;
    propensities[update.index] = update.propensity;
} // update()

/*---------------------------------------------------------------------------*/

void LinearSolver::update(std::vector<Update> updates)
{
    for (Update u : updates)
    {
        update(u);
    }
} // update()

/*---------------------------------------------------------------------------*/

std::optional<Event> LinearSolver::event()
{
    if (number_of_active_indices == 0)
    {
        propensity_sum = 0.0;
        return std::optional<Event>();
    }

    double r1 = sampler.generate();
    double r2 = sampler.generate();
    double fraction = propensity_sum * r1;
    double partial = 0.0;

    unsigned long m;

    for (m = 0; m < propensities.size(); m++)
    {
        partial += propensities[m];
        if (partial > fraction)
            break;
    }

    double dt = -std::log(r2) / propensity_sum;
    if (m < propensities.size())
        return std::optional<Event>(Event{.index = m, .dt = dt});
    else
        return std::optional<Event>(Event{.index = last_non_zero_event, .dt = dt});
} // event()

/*---------------------------------------------------------------------------*/

double LinearSolver::get_propensity(int index)
{
    return propensities[index];
} // get_propensity()

/*---------------------------------------------------------------------------*/

double LinearSolver::get_propensity_sum()
{
    return propensity_sum;
} // get_propensity_sum()
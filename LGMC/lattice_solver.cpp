/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://lzichi.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#include <iostream>
#include <vector>

#include "lattice_solver.h"

LatticeSolver::LatticeSolver(unsigned long int seed,
                             std::vector<double> &&initial_propensities) : propensity_sum(0.0),
                                                                           number_of_active_indices(0),
                                                                           // if this move isn't here, the semantics is that initial
                                                                           // propensities gets moved into a stack variable for the function
                                                                           // call and that stack variable is copied into the object.
                                                                           propensities(std::move(initial_propensities)),
                                                                           sampler(Sampler(seed))
{
    for (unsigned long i = 0; i < propensities.size(); i++)
    {
        propensity_sum += propensities[i];
        if (propensities[i] > 0)
        {
            number_of_active_indices += 1;
        }
    }
} // LatticeSolver()

/* ---------------------------------------------------------------------- */

LatticeSolver::LatticeSolver(unsigned long int seed,
                             std::vector<double> &initial_propensities) : propensity_sum(0.0),
                                                                          number_of_active_indices(0),
                                                                          propensities(initial_propensities),
                                                                          sampler(Sampler(seed))
{
    for (unsigned long i = 0; i < propensities.size(); i++)
    {
        propensity_sum += propensities[i];
        if (propensities[i] > 0)
        {
            number_of_active_indices += 1;
        }
    }
} // LatticeSolver()

/* ---------------------------------------------------------------------- */

void LatticeSolver::update(Update update)
{

    if (propensities[update.index] > 0.0)
        number_of_active_indices--;

    if (update.propensity > 0.0)
    {
        number_of_active_indices++;
    }

    propensity_sum -= propensities[update.index];
    propensity_sum += update.propensity;

    propensities[update.index] = update.propensity;
} // update()

/* ---------------------------------------------------------------------- */

void LatticeSolver::update(LatticeUpdate lattice_update, std::unordered_map<std::string,
                                                                            std::vector<std::pair<double, int>>> &props)
{

    propensity_sum += lattice_update.propensity;
    number_of_active_indices++;

    std::string hash = make_string(lattice_update.site_one, lattice_update.site_two);
    props[hash].push_back(std::make_pair(lattice_update.propensity, lattice_update.index));
} // update()

/* ---------------------------------------------------------------------- */

void LatticeSolver::update(std::vector<LatticeUpdate> lattice_updates,
                           std::unordered_map<std::string,
                                              std::vector<std::pair<double, int>>> &props)
{
    for (LatticeUpdate u : lattice_updates)
    {
        update(u, props);
    }
} // update()

/* ---------------------------------------------------------------------- */

void LatticeSolver::update(std::vector<Update> updates)
{
    for (Update u : updates)
    {
        update(u);
    }
} // update()

/* ---------------------------------------------------------------------- */

std::optional<LatticeEvent> LatticeSolver::event_lattice(std::unordered_map<std::string,
                                                                            std::vector<std::pair<double, int>>> &props)
{
    if (!(propensity_sum > 0))
    {
        double sum = 0;
        for (int i = 0; i < int(propensities.size()); i++)
        {
            sum += propensities[i];
        }
        for (auto it = props.begin(); it != props.end(); it++)
        {
            for (int i = 0; i < int(it->second.size()); i++)
            {
                sum += it->second[i].first;
            }
        }
        propensity_sum = sum;

        if (!(propensity_sum > 0))
        {
            return std::optional<LatticeEvent>();
        }
    }
    if (number_of_active_indices == 0)
    {
        propensity_sum = 0.0;
        return std::optional<LatticeEvent>();
    }

    bool isFound = false;
    unsigned long m;
    unsigned long int reaction_id;
    std::optional<int> site_one;
    std::optional<int> site_two;
    std::string hash;
    double dt;

    while (!isFound)
    {
        reaction_id = 0;
        site_one = std::optional<int>();
        site_two = std::optional<int>();
        dt = 0.;

        double r1 = sampler.generate();
        double r2 = sampler.generate();
        double fraction = propensity_sum * r1;
        long double partial = 0.0;

        // start with Gillespie propensities
        for (m = 0; m < propensities.size(); m++)
        {
            partial += propensities[m];
            if (partial > fraction)
            {
                isFound = true;
                reaction_id = m;
                break;
            }
        }

        // go through lattice propensities if not found
        if (!isFound)
        {
            auto it = props.begin();
            while (!isFound && it != props.end())
            {
                for (int i = 0; i < int(it->second.size()); i++)
                {

                    partial += it->second[i].first;

                    if (partial > fraction)
                    {

                        isFound = true;
                        hash = it->first;
                        reaction_id = it->second[i].second;

                        std::size_t pos = hash.find(".");
                        site_one = std::optional<int>(stoi(hash.substr(0, pos)));
                        site_two = std::optional<int>(stoi(hash.substr(pos + 1)));
                        if (site_one < site_two)
                        {
                            assert(false);
                        }
                        isFound = true;
                        break;
                    }
                }
                it++;
            }
        }

        dt = -std::log(r2) / propensity_sum;

        // check if found, if not propensity sum is incorrect
        if (!isFound || !(propensity_sum > 0))
        {
            double sum = 0;
            for (int i = 0; i < static_cast<int>(propensities.size()); i++)
            {
                sum += propensities[i];
            }
            for (auto it = props.begin(); it != props.end(); it++)
            {
                for (int i = 0; i < int(it->second.size()); i++)
                {
                    sum += it->second[i].first;
                }
            }
            propensity_sum = sum;

            if (!(propensity_sum > 0))
            {
                return std::optional<LatticeEvent>();
            }
        }
    }

    return std::optional<LatticeEvent>(LatticeEvent{.site_one = site_one, .site_two = site_two, .index = reaction_id, .dt = dt});
} // event_lattice()

/* ---------------------------------------------------------------------- */

std::string LatticeSolver::make_string(int site_one, int site_two)
{
    return (site_one > site_two) ? std::to_string(site_one) + "." +
                                       std::to_string(site_two)
                                 : std::to_string(site_two) + "." + std::to_string(site_one);

} // make_string()

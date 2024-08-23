/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#include "energy_reaction_network_simulation.h"

template <typename Solver>
void EnergyReactionNetworkSimulation<Solver>::init()
{
    std::vector<double> initial_propensities_temp;

    energy_reaction_network.compute_initial_propensities(state.homogeneous, initial_propensities_temp);
    solver = Solver(this->seed, std::ref(initial_propensities_temp));
    this->update_function = [&](Update update)
    { solver.update(update); };
} // init()

/* ------------------------------------------------------------------- */

template <typename Solver>
bool EnergyReactionNetworkSimulation<Solver>::execute_step()
{
    std::optional<Event> maybe_event = solver.event();

    if (!maybe_event)
    {

        return false;
    }
    else
    {
        // an event happens
        Event event = maybe_event.value();
        int next_reaction = event.index;

        // update time
        this->time += event.dt;

        // record what happened
        history.push_back(ReactionNetworkTrajectoryHistoryElement{
            .seed = this->seed,
            .reaction_id = next_reaction,
            .time = this->time,
            .step = this->step});

        if (history.size() == this->history_chunk_size)
        {
            history_queue.insert_history(
                std::move(
                    HistoryPacket<ReactionNetworkTrajectoryHistoryElement>{
                        .seed = this->seed,
                        .history = std::move(this->history)}));

            history = std::vector<ReactionNetworkTrajectoryHistoryElement>();
            history.reserve(this->history_chunk_size);
        }

        // increment step
        this->step++;

        // update state
        energy_reaction_network.update_state(std::ref(state.homogeneous), next_reaction);

        // update propensities
        energy_reaction_network.update_energy_budget(std::ref(state.energy_budget), next_reaction);

        // update propensities
        energy_reaction_network.update_propensities(
            this->update_function,
            std::ref(state.homogeneous),
            next_reaction,
            state.energy_budget);

        return true;
    }
} // execute_step()

/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://lzichi.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#include "reaction_network_simulation.h"
#include "../GMC/reaction_network.h"

template <typename Solver>
void ReactionNetworkSimulation<Solver>::init()
{
    std::vector<double> initial_propensities_temp;
    reaction_network.compute_initial_propensities(state, initial_propensities_temp);
    solver = Solver(this->seed, std::ref(initial_propensities_temp));
    this->update_function = [&](Update update)
    { solver.update(update); };

} // init()

/* ------------------------------------------------------------------- */

template <typename Solver>
bool ReactionNetworkSimulation<Solver>::execute_step()
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
        reaction_network.update_state(std::ref(state), next_reaction);

        // update propensities
        reaction_network.update_propensities(
            this->update_function,
            std::ref(state),
            next_reaction);

        return true;
    }
} // execute_step()

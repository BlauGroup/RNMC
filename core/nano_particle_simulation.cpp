/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://lzichi.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#include "nano_particle_simulation.h"

void NanoParticleSimulation::init()
{
    std::vector<std::set<int>> seed_site_reaction_dependency;
    std::vector<NanoReaction> seed_reactions;

    seed_site_reaction_dependency.resize(nano_particle.sites.size());
    nano_particle.compute_reactions(state, std::ref(seed_reactions), std::ref(seed_site_reaction_dependency));
    nanoSolver = NanoSolver(this->seed, std::ref(seed_reactions));
    site_reaction_dependency = seed_site_reaction_dependency;
} // init()

/* ------------------------------------------------------------------- */

bool NanoParticleSimulation::execute_step()
{
    std::optional<Event> maybe_event = nanoSolver.event();

    if (!maybe_event)
    {

        return false;
    }
    else
    {
        // an event happens
        Event event = maybe_event.value();
        int next_reaction_id = event.index;
        NanoReaction next_reaction = nanoSolver.current_reactions[next_reaction_id];

        // update time
        this->time += event.dt;

        // record what happened
        history.push_back(NanoTrajectoryHistoryElement{
            .seed = this->seed,
            .reaction = next_reaction,
            .time = this->time,
            .step = this->step});

        if (history.size() == this->history_chunk_size)
        {
            history_queue.insert_history(
                HistoryPacket<NanoTrajectoryHistoryElement>{
                    .seed = this->seed,
                    .history = std::move(this->history)});

            history = std::vector<NanoTrajectoryHistoryElement>();
            history.reserve(this->history_chunk_size);
        }

        // increment step
        this->step++;

        // update state
        nano_particle.update_state(std::ref(state), next_reaction);

        // update list of current available reactions
        nano_particle.update_reactions(std::cref(state), next_reaction,
                                       std::ref(site_reaction_dependency),
                                       std::ref(nanoSolver.current_reactions));
        nanoSolver.update();

        return true;
    }
} // execute_step()
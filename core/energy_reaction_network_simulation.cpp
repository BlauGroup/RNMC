#include "energy_reaction_network_simulation.h"
#include "../GMC/energy_reaction_network.h"

template <typename Solver>
void EnergyReactionNetworkSimulation<Solver>::init() {
    energy_reaction_network.compute_initial_propensities(state);
    solver = Solver(this->seed, std::ref(energy_reaction_network.initial_propensities));
    this->update_function = [&] (Update update) {solver.update(update);};
    energy_budget = energy_reaction_network.energy_budget;
} // init()

/* ---------------------------------------------------------------------- */

template <typename Solver>
bool EnergyReactionNetworkSimulation<Solver>::execute_step() {
    std::optional<Event> maybe_event = solver.event();

    if (! maybe_event) {

        return false;

    } else {
        // an event happens
        Event event = maybe_event.value();
        int next_reaction = event.index;

        // update time
        this->time += event.dt;

        // record what happened
        history.push_back(ReactionNetworkTrajectoryHistoryElement {
            .seed = this->seed,
            .reaction_id = next_reaction,
            .time = this->time,
            .step = this->step
            });

        if (history.size() == this->history_chunk_size ) {
            history_queue.insert_history(
                std::move(
                    HistoryPacket<ReactionNetworkTrajectoryHistoryElement> {
                        .history = std::move(this->history),
                        .seed = this->seed
                    }));

            history = std::vector<ReactionNetworkTrajectoryHistoryElement> ();
            history.reserve(this->history_chunk_size);
        }

        // increment step
        this->step++;

        // update state
        energy_reaction_network.update_state(std::ref(state), next_reaction);

        // update propensities
        energy_reaction_network.update_energy_budget(std::ref(energy_budget), next_reaction);
            
        // update propensities
        energy_reaction_network.update_propensities(
            update_function,
            std::ref(state),
            next_reaction,
            energy_budget);

        return true;
    }
} // execute_step()

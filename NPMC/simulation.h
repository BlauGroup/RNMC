#pragma once
#include "nano_particle.h"
#include "../core/solvers.h"


struct HistoryElement {
    int site_id_1;
    int site_id_2;
    int interaction_id;
    double time;
};

template <typename Solver>
struct Simulation {
    NanoParticle &nano_particle;
    unsigned long int seed;
    std::vector<int> state;
    double time;
    int step;       // number of reactions which have occoured
    Solver solver;
    std::vector<HistoryElement> history;


    Simulation( NanoParticle &nano_particle,
                unsigned long int seed,
                int step_cutoff
        ) :
        nano_particle (nano_particle),
        seed (seed),
        state (nano_particle.initial_state),
        time (0.0),
        step (0),
        solver (seed, std::ref(nano_particle.initial_propensities)),
        history (step_cutoff + 1)
        {};

    bool execute_step();
    void execute_steps(int step_cutoff);
};

template <typename Solver>
bool Simulation<Solver>::execute_step() {

    std::optional<Event> maybe_event = solver.event();
    if (! maybe_event) {
        return false;
    } else {
        // an event happens
        Event event = maybe_event.value();
        int next_reaction_id = event.index;
        Reaction reaction = nano_particle.reactions[next_reaction_id];

        // update time
        time += event.dt;

        history[step] = HistoryElement {
            .site_id_1 = reaction.site_id_1,
            .site_id_2 = reaction.site_id_2,
            .interaction_id = reaction.interaction_id,
            .time = time};

        // increment step
        step++;


        // update state
        Interaction interaction = nano_particle.interactions[
            reaction.interaction_id];

        state[reaction.site_id_1] = interaction.right_state_1;
        if (interaction.number_of_sites == 2) {
            state[reaction.site_id_2] = interaction.right_state_2;
        }
        // update propensities

        for ( unsigned int i = 0;
              i < nano_particle.site_reaction_dependency[reaction.site_id_1].size();
              i++ ) {
            int reaction_id = nano_particle.site_reaction_dependency[reaction.site_id_1][i];
            double new_propensity = nano_particle.compute_propensity(std::ref(state), reaction_id);


            solver.update( Update {
                    .index = (unsigned long int) reaction_id,
                    .propensity = new_propensity});
        }

        if (interaction.number_of_sites == 2) {
            for ( unsigned int i = 0;
                  i < nano_particle.site_reaction_dependency[reaction.site_id_2].size();
                  i++ ) {
                int reaction_id = nano_particle.site_reaction_dependency[reaction.site_id_2][i];
                double new_propensity = nano_particle.compute_propensity(std::ref(state), reaction_id);

                solver.update( Update {
                        .index = (unsigned long int) reaction_id,
                        .propensity = new_propensity});
            }

        }
        return true;

    }
}

template <typename Solver>
void Simulation<Solver>::execute_steps(int step_cutoff) {
    while(execute_step()) {
        if (step > step_cutoff)
            break;
    }
};

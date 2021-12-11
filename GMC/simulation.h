#pragma once
#include "reaction_network.h"
#include "../core/solvers.h"
#include <functional>

struct HistoryElement {
    int reaction; // reaction which fired
    double time;  // time after reaction has occoured.
};

template <typename Solver>
struct Simulation {
    ReactionNetwork &reaction_network;
    unsigned long int seed;
    std::vector<int> state;
    double time;
    int step; // number of reactions which have occoured
    Solver solver;
    std::vector<HistoryElement> history;


    Simulation(ReactionNetwork &reaction_network,
               unsigned long int seed,

               // step cutoff gets used here to set the history length
               // we don't actually store it in the Simulation object
               int step_cutoff) :
        reaction_network (reaction_network),
        seed (seed),
        state (reaction_network.initial_state),
        time (0.0),
        step (0),
        solver (seed, std::ref(reaction_network.initial_propensities)),
        history (step_cutoff + 1)
        {};


    bool execute_step();
    void execute_steps(int step_cutoff);
    bool check_state_positivity();
};

template <typename Solver>
bool Simulation<Solver>::check_state_positivity() {
    for (int i = 0; i < state.size(); i++) {
    if (state[i] < 0) {
        std::cerr << time_stamp()
                  << "negative state encountered: index = " << i << "\n";
      return false;
    }
  }
  return true;
};


template <typename Solver>
bool Simulation<Solver>::execute_step() {
    std::optional<Event> maybe_event = solver.event();

    if (! maybe_event) {

        return false;

    } else {
        // an event happens
        Event event = maybe_event.value();
        int next_reaction = event.index;

        // update time
        time += event.dt;

        // record what happened
        history[step] = HistoryElement {
            .reaction = next_reaction,
            .time = time};

        // increment step
        step++;

        // update state
        for (int m = 0;
             m < reaction_network.reactions[next_reaction].number_of_reactants;
             m++) {
            state[reaction_network.reactions[next_reaction].reactants[m]]--;
        }

        for (int m = 0;
             m < reaction_network.reactions[next_reaction].number_of_products;
             m++) {
            state[reaction_network.reactions[next_reaction].products[m]]++;
        }

        // update propensities
        std::optional<std::vector<int>> &maybe_dependents =
            reaction_network.get_dependency_node(next_reaction);

        if (maybe_dependents) {
            // relevent section of dependency graph has been computed
            std::vector<int> &dependents = maybe_dependents.value();

            for (unsigned long int m = 0; m < dependents.size(); m++) {
                unsigned long int reaction_index = dependents[m];
                double new_propensity = reaction_network.compute_propensity(
                    state,
                    reaction_index);

                solver.update(Update {
                        .index = reaction_index,
                        .propensity = new_propensity});

            }
        } else {
            // relevent section of dependency graph has not been computed
            for (unsigned long int reaction_index = 0;
                 reaction_index < reaction_network.reactions.size();
                 reaction_index++) {

                double new_propensity = reaction_network.compute_propensity(
                    state,
                    reaction_index);

                solver.update(Update {
                        .index = reaction_index,
                        .propensity = new_propensity});

            }

        }

        return true;
    }
};

template <typename Solver>
void Simulation<Solver>::execute_steps(int step_cutoff) {
    while(execute_step()) {
        // check_state_positivity();
        if (step > step_cutoff)
            break;
    }
};

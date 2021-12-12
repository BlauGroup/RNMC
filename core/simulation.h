#pragma once
#include "solvers.h"
#include <functional>


struct HistoryElement {
    int reaction_id; // reaction which fired
    double time;  // time after reaction has occoured.
};

template <typename Solver, typename Model>
struct Simulation {
    Model &model;
    unsigned long int seed;
    std::vector<int> state;
    double time;
    int step; // number of reactions which have occoured
    Solver solver;
    std::vector<HistoryElement> history;


    Simulation(Model &model,
               unsigned long int seed,

               // step cutoff gets used here to set the history length
               // we don't actually store it in the Simulation object
               int step_cutoff) :
        model (model),
        seed (seed),
        state (model.initial_state),
        time (0.0),
        step (0),
        solver (seed, std::ref(model.initial_propensities)),
        history (step_cutoff + 1)
        {};


    bool execute_step();
    void execute_steps(int step_cutoff);
};


template <typename Solver, typename Model>
bool Simulation<Solver, Model>::execute_step() {
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
            .reaction_id = next_reaction,
            .time = time};

        // increment step
        step++;

        // update state
        model.update_state(std::ref(state), next_reaction);


        // update propensities
        model.update_propensities(
            std::ref(solver),
            std::ref(state),
            next_reaction);

        return true;
    }
};

template <typename Solver, typename Model>
void Simulation<Solver, Model>::execute_steps(int step_cutoff) {
    while(execute_step()) {
        if (step > step_cutoff)
            break;
    }
};

#pragma once
#include "solvers.h"
#include "queues.h"
#include <functional>
#include <csignal>

struct HistoryElement {

    unsigned long int seed; // seed
    int reaction_id; // reaction which fired
    double time;  // time after reaction has occoured.
    int step;
};

struct HistoryPacket {
    std::vector<HistoryElement> history;
    unsigned long int seed;
};

template <typename Solver, typename Model>
struct Simulation {
    Model &model;
    unsigned long int seed;
    std::vector<int> state;
    double time;
    int step; // number of reactions which have occoured
    Solver solver;
    unsigned long int history_length;
    std::vector<HistoryElement> history;
    HistoryQueue<HistoryPacket> &history_queue;
    std::function<void(Update)> update_function;


    Simulation(Model &model,
               unsigned long int seed,
               int history_length,
               HistoryQueue<HistoryPacket> &history_queue
        ) :
        model (model),
        seed (seed),
        state (model.initial_state),
        time (0.0),
        step (0),
        solver (seed, std::ref(model.initial_propensities)),
        history_length(history_length),
        history_queue(history_queue),
        update_function ([&] (Update update) {solver.update(update);})
        {
            history.reserve(history_length);
        };


    bool execute_step();
    void execute_steps(int step_cutoff);
    void execute_time(double time_cutoff);

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
        history.push_back(HistoryElement {
            .seed = seed,
            .reaction_id = next_reaction,
            .time = time,
            .step = step
            });

        if ( history.size() == history_length ) {
            raise(SIGINT);
            history_queue.insert_history(
                std::move(
                    HistoryPacket {
                        .history = std::move(history),
                        .seed = seed
                        }));

            history = std::vector<HistoryElement> ();
            history.reserve(history_length);
        }


        // increment step
        step++;

        // update state
        model.update_state(std::ref(state), next_reaction);


        // update propensities
        model.update_propensities(
            update_function,
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

template <typename Solver, typename Model>
void Simulation<Solver, Model>::execute_time(double time_cutoff) {
    while(execute_step()) {
        if (time > time_cutoff)
            break;
    }
};

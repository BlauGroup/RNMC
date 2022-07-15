#pragma once
#include "solvers.h"
#include "../LGMC/LatSolver.h"
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
    std::unordered_map<std::string,                     
        std::vector< std::pair<double, int> > > props;
    double time;
    int step; // number of reactions which have occoured
    Solver solver;
    unsigned long int history_chunk_size;
    std::vector<HistoryElement> history;
    HistoryQueue<HistoryPacket> &history_queue;
    std::function<void(Update)> update_function;
    std::function<void(LatticeUpdate)> lattice_update_function;

    Simulation(Model &model,
               unsigned long int seed,
               int history_chunk_size,
               HistoryQueue<HistoryPacket> &history_queue
        ) :
        model (model),
        seed (seed),
        state (model.initial_state),
        time (0.0),
        step (0),
        solver (seed, std::ref(model.initial_propensities)),
        history_chunk_size (history_chunk_size),
        history_queue(history_queue),
        update_function ([&] (Update update) {solver.update(update);})
        {
            history.reserve(history_chunk_size);
        };


    bool execute_step();
    bool execute_step_LGMC();
    void init_lattice_props(); // TODO: FOR EVAN
    void execute_steps(int step_cutoff, bool isLGMC);
    void execute_time(double time_cutoff, bool isLGMC);

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

        if ( history.size() == history_chunk_size ) {
            history_queue.insert_history(
                std::move(
                    HistoryPacket {
                        .history = std::move(history),
                        .seed = seed
                        }));

            history = std::vector<HistoryElement> ();
            history.reserve(history_chunk_size);
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
bool Simulation<Solver, Model>::execute_step_LGMC() {
    std::pair<std::optional<Event>, std::optional<LatticeEvent>> maybe_events = solver.event_lattice();

    int next_reaction = 0;

    if (!maybe_events.first && !maybe_events.second) {

        return false;

    } 
    else if (maybe_events.second) {
        // lattice event happens
        LatticeEvent event = maybe_events.second.value();
        int next_reaction = event.index;

        // update time
        time += event.dt;

        // update state
        bool update_gillepsie = model.update_state(std::ref(props), next_reaction, 
                                                    event.site_one, event.site_two, solver.propensity_sum);
        model.update_propensities(std::ref(state), lattice_update_function, next_reaction, event.site_one, event.site_two);
        
        if(update_gillepsie) {
            model.update_state(std::ref(state), next_reaction);
            model.update_propensities(update_function, std::ref(state), next_reaction);
        }
    }
    
    else {
        // gillespie event happens
        Event event = maybe_events.first.value();
        int next_reaction = event.index;

        // update time
        time += event.dt;

        // update state
        model.update_state(std::ref(state), next_reaction);

        // update propensities
        model.update_propensities(
            update_function,
            std::ref(state),
            next_reaction);
    }

    // record what happened
    history.push_back(HistoryElement {
        .seed = seed,
        .reaction_id = next_reaction,
        .time = time,
        .step = step
        });

    if ( history.size() == history_chunk_size ) {
        history_queue.insert_history(
            std::move(
                HistoryPacket {
                    .history = std::move(history),
                    .seed = seed
                    }));

        history = std::vector<HistoryElement> ();
        history.reserve(history_chunk_size);
    }


    // increment step
    step++;

    return true;
};

template <typename Solver, typename Model>
void Simulation<Solver, Model>::execute_steps(int step_cutoff, bool isLGMC) {
    
    
    init_lattice_props(); // TODO: FOR EVAN

    if(isLGMC) {
        while(execute_step_LGMC()) {
        if (step > step_cutoff)
            break;
        }
    }
    else {
        while(execute_step()) {
            if (step > step_cutoff)
                break;
        }
    }
    
};

template <typename Solver, typename Model>
void Simulation<Solver, Model>::execute_time(double time_cutoff, bool isLGMC) {
    if(isLGMC) {
        while(execute_step_LGMC()) {
        if (time > time_cutoff)
            break;
        }
    }
    else {
        while(execute_step()) {
        if (time > time_cutoff)
            break;
        }
    }
};

/* ----------------*/
/* TODO: FOR EVAN */
/* ----------------*/
template <typename Solver, typename Model>
void Simulation<Solver, Model>::init_lattice_props() {
    assert(false);
}
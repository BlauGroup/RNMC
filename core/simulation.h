#pragma once
#include "solvers.h"
#include "queues.h"
#include <functional>
#include <csignal>
#include <set>
//
// //Include this here for now, but breakout into own file later
// struct Reaction {
//     int site_id[2];
//     Interaction interaction;
//
//     // rate has units 1 / s
//     double rate;
// };

namespace {
  // In the GNUC Library, sig_atomic_t is a typedef for int,
  // which is atomic on all systems that are supported by the
  // GNUC Library
  volatile sig_atomic_t do_shutdown = 0;

  // std::atomic is safe, as long as it is lock-free
  std::atomic<bool> shutdown_requested = false;
  static_assert( std::atomic<bool>::is_always_lock_free );
  // or, at runtime: assert( shutdown_requested.is_lock_free() );
}

struct CutoffHistoryElement{
    unsigned long int seed;
    int step;
    double time;
};

struct CutoffHistoryPacket{
    unsigned long int seed;
    std::vector<CutoffHistoryElement> cutoff;
};

struct StateHistoryElement{
    unsigned long int seed; //seed
    int site_id; 
    int degree_of_freedom; //energy level the site is at
};

struct StateHistoryPacket {
    unsigned long int seed; //seed
    std::vector<StateHistoryElement> state; // current state of the reaction to be written
};

struct HistoryElement {
    unsigned long int seed; // seed
    Reaction reaction; // reaction which fired
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
    unsigned long int history_chunk_size;
    std::vector<HistoryElement> history;
    HistoryQueue<HistoryPacket> &history_queue;
    HistoryQueue<StateHistoryPacket> &state_history_queue;
    std::vector<std::set<int>> site_reaction_dependency;


    Simulation(Model &model,
               unsigned long int seed,
               int step,
               double time,
               std::vector<int> state,
               std::vector<Reaction> initial_reactions,
               std::vector<std::set<int>> site_reaction_dependency,
               int history_chunk_size,
               HistoryQueue<HistoryPacket> &history_queue,
               HistoryQueue<StateHistoryPacket> &state_history_queue
        ) :
        model (model),
        seed (seed),
        state (state),
        time (time),
        step (step),
        solver (seed, std::ref(initial_reactions)),
        history_chunk_size (history_chunk_size),
        history_queue (history_queue),
        state_history_queue (state_history_queue),
        // update_function ([&] (Update update) {solver.update(update);}),
        // current_reactions (model.initial_reactions),
        site_reaction_dependency (site_reaction_dependency)
        {
            history.reserve(history_chunk_size);
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
        int next_reaction_id = event.index;
        Reaction next_reaction = solver.current_reactions[next_reaction_id];

        // update time
        time += event.dt;

        // record what happened
        history.push_back(HistoryElement {
            .seed = seed,
            .reaction = next_reaction,
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

        // update list of current available reactions
        model.update_reactions(std::cref(state), next_reaction, std::ref(site_reaction_dependency), std::ref(solver.current_reactions));
        solver.update();

        return true;
    }
};

void write_error_message(std::string s){
    char char_array[s.length()+1];
    strcpy(char_array, s.c_str());

    write(STDERR_FILENO, char_array, sizeof(char_array) - 1);
}

template <typename Solver, typename Model>
void Simulation<Solver, Model>::execute_steps(int step_cutoff) {
    while(execute_step()) {
        if (step > step_cutoff) {
            break;
        } else if (do_shutdown || shutdown_requested.load()) {
            // Handle shutdown request from SIGTERM
            write_error_message("Received termination request on thread - cleaning up\n");
            break;
        }
    }
};

template <typename Solver, typename Model>
void Simulation<Solver, Model>::execute_time(double time_cutoff) {
    while(execute_step()) {
        if (time > time_cutoff) {
            break;
        } else if (do_shutdown || shutdown_requested.load()) {
            // Handle shutdown request from SIGTERM
            write_error_message("Received termination request on thread - cleaning up\n");
            break;
        }
    }
};

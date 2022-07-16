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
class Simulation {
    private: 
        Solver solver;
    public:
    Model &model;
    unsigned long int seed;
    std::vector<int> state;
    double time;
    int step; // number of reactions which have occoured
    unsigned long int history_chunk_size;
    std::vector<HistoryElement> history;
    HistoryQueue<HistoryPacket> &history_queue;
    std::function<void(Update)> update_function;

    Simulation(Model &model, unsigned long int seed,
               int history_chunk_size,
               HistoryQueue<HistoryPacket> &history_queue
        ) :
        model (model),
        seed (seed),
        state (model.initial_state),
        time (0.0),
        step (0),
        history_chunk_size (history_chunk_size),
        history_queue(history_queue),
        update_function ([&] (Update update) {solver.update(update);})
        {
            history.reserve(history_chunk_size);
        };

    void init(Model &model);
    virtual bool execute_step();
    void execute_steps(int step_cutoff);
    void execute_time(double time_cutoff);

};
template <typename Solver, typename Model>
void Simulation<Solver, Model>::init(Model &model) {
    solver = Solver(seed, std::ref(model.initial_propensities));
}

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

/* ---------------------------------------------------------------------- */
template<typename Model>
class LatticeSimulation : public Simulation<LatSolver, Model> {
    public:
    std::unordered_map<std::string,                     
        std::vector< std::pair<double, int> > > props;
    LatSolver latsolver;
    std::function<void(LatticeUpdate)> lattice_update_function;

    LatticeSimulation(Model &model, unsigned long int seed, int history_chunk_size,
               HistoryQueue<HistoryPacket> &history_queue) :
               Simulation<LatSolver, Model>(model, seed, history_chunk_size, history_queue),
               latsolver (seed, std::ref(model.initial_propensities)),
               lattice_update_function ([&] (LatticeUpdate lattice_update) {latsolver.update(lattice_update);}) {};

    bool execute_step();
    void init(); 

};

template<typename Model>
void LatticeSimulation<Model>::init() {
    // TODO: initialize lattice components
    assert(false);
}

template<typename Model>
bool LatticeSimulation<Model>::execute_step() {

    std::optional<LatticeEvent> maybe_event = latsolver.event_lattice();

    if (! maybe_event) {

        return false;

    } else {
        // an event happens
        LatticeEvent event = maybe_event.value();
        int next_reaction = event.index;

        // update time
        this->time += event.dt;

        // record what happened
        this->history.push_back(HistoryElement {
            .seed = this->seed,
            .reaction_id = next_reaction,
            .time = this->time,
            .step = this->step
            });

        if (this->history.size() == this->history_chunk_size ) {
            this->history_queue.insert_history(
                std::move(
                    HistoryPacket {
                        .history = std::move(this->history),
                        .seed = this->seed
                        }));

            this->history = std::vector<HistoryElement> ();
            this->history.reserve(this->history_chunk_size);
        }


        // increment step
        this->step++;

        // update_state
        this->model.update_state(std::ref(props), std::ref(this->state), next_reaction, 
                    event.site_one, event.site_two, latsolver.propensity_sum);

        // update_propensities 
        this->model.update_propensities(std::ref(this->state), this->update_function, 
                                        lattice_update_function, next_reaction, 
                                        event.site_one, event.site_two);

        return true;
    }
 
}

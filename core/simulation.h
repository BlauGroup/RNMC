#pragma once
#include "solvers.h"
#include "../LGMC/LatSolver.h"
#include "queues.h"
#include "../LGMC/lattice.h"
#include <functional>
#include <csignal>

struct HistoryElement {

    unsigned long int seed; // seed
    int reaction_id; // reaction which fired
    double time;  // time after reaction has occoured.
    int step;
};

template <typename T>
struct HistoryPacket {
    std::vector<T> history;
    unsigned long int seed;
};

struct LatticeHistoryElement {
    unsigned long int seed;
    int step;
    int reaction_id;
    int site_1;
    int site_2;
    double time;
};

template <typename Solver, typename Model, typename History>
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
    std::vector<History> history;
    HistoryQueue<HistoryPacket<History>> &history_queue;   
    std::function<void(Update)> update_function;


    Simulation(Model &model, unsigned long int seed,
               int history_chunk_size,
               HistoryQueue<HistoryPacket<History>> &history_queue
        ) :
        model (model),
        seed (seed),
        state (model.initial_state),
        history_queue(history_queue),
        time (0.0),
        step (0),
        history_chunk_size (history_chunk_size)
        {
            history.reserve(history_chunk_size);
        };

    void init(Model &model);
    virtual bool execute_step();
    void execute_steps(int step_cutoff);
    void execute_time(double time_cutoff);

};
template <typename Solver, typename Model, typename History>
void Simulation<Solver, Model, History>::init(Model &model) {
    solver = Solver(seed, std::ref(model.initial_propensities));
    update_function = std::function<void(Update)> ([&] (Update update) {solver.update(update);});
    history.reserve(history_chunk_size);
}

template <typename Solver, typename Model, typename History>
bool Simulation<Solver, Model, History>::execute_step() {
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
        history.push_back(History {
            .seed = seed,
            .reaction_id = next_reaction,
            .time = time,
            .step = step
            });

        if ( history.size() == history_chunk_size ) {
            history_queue.insert_history(
                std::move(
                    HistoryPacket<History> {
                        .history = std::move(history),
                        .seed = seed
                        }));

            history = std::vector<History> ();
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

template <typename Solver, typename Model, typename History>
void Simulation<Solver, Model, History>::execute_steps(int step_cutoff) {

    while(execute_step()) {
        if (step > step_cutoff)
            break;
    }
    
};

template <typename Solver, typename Model, typename History>
void Simulation<Solver, Model, History>::execute_time(double time_cutoff) {
    while(execute_step()) {
        if (time > time_cutoff)
            break;
    }
};

/* ---------------------------------------------------------------------- */
template<typename Model, typename History>
class LatticeSimulation : public Simulation<LatSolver, Model, History> {
    public:
    std::unordered_map<std::string,                     
        std::vector< std::pair<double, int> > > props;
    LatSolver latsolver;
    Lattice *lattice;
    std::function<void(LatticeUpdate, std::unordered_map<std::string,                     
                std::vector< std::pair<double, int> > > &)> lattice_update_function;
    std::function<void(Update)> update_function;

    LatticeSimulation(Model &model, unsigned long int seed, int history_chunk_size,
               HistoryQueue<HistoryPacket<History>> &history_queue) :
               Simulation<LatSolver, Model, History>(model, seed, history_chunk_size, history_queue),
               latsolver (seed, std::ref(model.initial_propensities)),
               lattice_update_function ([&] (LatticeUpdate lattice_update, 
               std::unordered_map<std::string,                     
                std::vector< std::pair<double, int> > > &props) {latsolver.update(lattice_update, props);}), 
               update_function ([&] (Update update) {latsolver.update(update);}) {
                
                    if(model.initial_lattice) {
                        lattice = new  Lattice(*model.initial_lattice);
                    } else {
                        lattice = nullptr;
                    }
                
               };

    bool execute_step();
    void init(); 

};

template<typename Model, typename History>
void LatticeSimulation<Model, History>::init() {
    this->model.update_adsorp_state(this->lattice, this->props, latsolver.propensity_sum, latsolver.number_of_active_indices);
    this->model.update_adsorp_props(this->lattice, lattice_update_function, this->state, std::ref(props));
    
}

template<typename Model, typename History>
bool LatticeSimulation<Model, History>::execute_step() {

    std::optional<LatticeEvent> maybe_event = latsolver.event_lattice(props);

    if (! maybe_event) {

        return false;

    } else {
        // an event happens
        LatticeEvent event = maybe_event.value();
        int next_reaction = event.index;

        // update time
        this->time += event.dt;
        int site_1 = -1;
        int site_2 = -1;

        if(event.site_one) {
            site_1 = event.site_one.value();
        }
        if(event.site_two) {
            site_2 = event.site_two.value();
        }

        if(site_1 == 5758) {
            std::cout << "next reaction: " << next_reaction << ", step: " << this->step << ", site 1: " << site_1 << ", site_2: " << site_2 << std::endl;
        }
        // record what happened
        this->history.push_back(History {
            .seed = this->seed,
            .reaction_id = next_reaction,
            .time = this->time,
            .step = this->step,
            .site_1 = site_1,
            .site_2 = site_2
            });
        if(this->history.size() == 1000) {
            std::cout << "1000 events\n" ;
        }

        if (this->history.size() == this->history_chunk_size ) {
            this->history_queue.insert_history(
                std::move(
                    HistoryPacket<History> {
                        .history = std::move(this->history),
                        .seed = this->seed
                        }));

            this->history = std::vector<History> ();
            this->history.reserve(this->history_chunk_size);
        }


        // increment step
        this->step++;

        // update_state
        this->model.update_state(lattice, std::ref(props), std::ref(this->state), next_reaction, 
                    event.site_one, event.site_two, latsolver.propensity_sum, latsolver.number_of_active_indices);

        // update_propensities 
        this->model.update_propensities(lattice, std::ref(this->state), this->update_function, 
                                        lattice_update_function, next_reaction, 
                                        event.site_one, event.site_two, props);
        return true;
    }
 
}

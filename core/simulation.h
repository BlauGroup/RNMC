#pragma once
#include "solvers.h"
#include "../LGMC/LatSolver.h"
#include "queues.h"
#include "../LGMC/lattice.h"
#include <functional>
#include <csignal>
#include <set>
#include <atomic>
#include <unistd.h>
#include <cstring>
#include "../GMC/reaction_network.h"
#include "../NPMC/NanoSolver.h"
#include "../NPMC/nano_particle.h"
#include "../LGMC/lattice_reaction_network.h"

void write_error_message(std::string s){
    char char_array[s.length()+1];
    strcpy(char_array, s.c_str());

    write(STDERR_FILENO, char_array, sizeof(char_array) - 1);
}

template <typename T>
struct HistoryPacket {
    std::vector<T> history;
    unsigned long int seed;
};


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

/*---------------------------------------------------------------------------*/

template <typename Solver>
class Simulation {
    public:
    unsigned long int seed;
    double time;
    int step; // number of reactions which have occoured
    unsigned long int history_chunk_size; 
    std::function<void(Update)> update_function;


    Simulation(unsigned long int seed,
               int history_chunk_size,
               int step,
               double time
        ) :
        seed (seed),
        time (time),
        step (step),
        history_chunk_size (history_chunk_size)
        {
        };

    void execute_steps(int step_cutoff);
    void execute_time(double time_cutoff);
    virtual bool execute_step();

};
template <typename Solver>
bool Simulation<Solver>::execute_step() {
    return false;
}

template <typename Solver>
void Simulation<Solver>::execute_steps(int step_cutoff) {

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

template <typename Solver>
void Simulation<Solver>::execute_time(double time_cutoff) {
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


/* ---------------------------------------------------------------------- */
template <typename Solver>
class ReactionNetworkSimulation : public Simulation<Solver> {
    private: 
        Solver solver;
    public:
    ReactionNetwork &reaction_network;
    std::vector<int> state;
    std::vector<ReactionNetworkTrajectoryHistoryElement> history;
    HistoryQueue<HistoryPacket<ReactionNetworkTrajectoryHistoryElement>> &history_queue; 

    ReactionNetworkSimulation(ReactionNetwork &reaction_network, 
            unsigned long int seed,
            int step,
            double time,
            std::vector<int> state,
            int history_chunk_size,
            HistoryQueue<HistoryPacket<ReactionNetworkTrajectoryHistoryElement>> &history_queue
        ) : 
        Simulation<Solver>(seed, history_chunk_size, step, time),
        history_queue(history_queue),
        reaction_network (reaction_network),
        state (state)
        { 
            history.reserve(this->history_chunk_size);
        };

    void init();
    bool execute_step();

};

template <typename Solver>
void ReactionNetworkSimulation<Solver>::init() {
    reaction_network.compute_initial_propensities();
    solver = Solver(this->seed, std::ref(reaction_network.initial_propensities));
    this->update_function = [&] (Update update) {solver.update(update);};

}

template <typename Solver>
bool ReactionNetworkSimulation<Solver>::execute_step() {
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
        reaction_network.update_state(std::ref(state), next_reaction);

        // update propensities
        reaction_network.update_propensities(
            this->update_function,
            std::ref(state),
            next_reaction);

        return true;
    }
};

/* ---------------------------------------------------------------------- */

class LatticeSimulation : public Simulation<LatSolver> {
    public:
    std::unordered_map<std::string,                     
        std::vector< std::pair<double, int> > > props;
    LatSolver latSolver;
    LatticeReactionNetwork &lattice_network;
    LatticeState state;


    std::function<void(LatticeUpdate, std::unordered_map<std::string,                     
                std::vector< std::pair<double, int> > > &)> lattice_update_function;
    std::function<void(Update)> update_function;

    std::vector<LatticeTrajectoryHistoryElement> history;
    HistoryQueue<HistoryPacket<LatticeTrajectoryHistoryElement>> &history_queue; 


    LatticeSimulation(LatticeReactionNetwork &lattice_network, unsigned long int seed, int step,
               double time, LatticeState state_in, int history_chunk_size,
               HistoryQueue<HistoryPacket<LatticeTrajectoryHistoryElement>> &history_queue) :
               Simulation<LatSolver>(seed, history_chunk_size, step, time),
               lattice_network(lattice_network),
               history_queue(history_queue)
                { 
                    state.homogeneous = state_in.homogeneous;
                    state.lattice = state_in.lattice;
                    state_in.lattice = NULL;
                    history.reserve(this->history_chunk_size);
     
                };

    bool execute_step();
    void init(); 
    ~LatticeSimulation() { delete state.lattice;};

};

void LatticeSimulation::init() {
    lattice_network.compute_initial_propensities(this->state.lattice);
    lattice_network.update_adsorp_state(this->state.lattice, this->props, latSolver.propensity_sum, latSolver.number_of_active_indices);
    lattice_network.update_adsorp_props(this->state.lattice, lattice_update_function, this->state.homogeneous, std::ref(props));

    latSolver = LatSolver(seed, std::ref(lattice_network.initial_propensities));
    this->update_function = [&] (Update update) {latSolver.update(update);};
    lattice_update_function = [&] (LatticeUpdate lattice_update, 
               std::unordered_map<std::string,                     
                std::vector< std::pair<double, int> > > &props) 
                 {latSolver.update(lattice_update, props);};

               
    
}

bool LatticeSimulation::execute_step() {

    std::optional<LatticeEvent> maybe_event = latSolver.event_lattice(props);

    if (!maybe_event) {

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

        if(history.size() == 1000) {
            std::cout << "1000 events" << std::endl;
        }
        // record what happened
        history.push_back(LatticeTrajectoryHistoryElement {
            .seed = this->seed,
            .reaction_id = next_reaction,
            .time = this->time,
            .step = this->step,
            .site_1 = site_1,
            .site_2 = site_2
            });

        if (history.size() == this->history_chunk_size ) {
            history_queue.insert_history(
                std::move(
                    HistoryPacket<LatticeTrajectoryHistoryElement> {
                        .history = std::move(this->history),
                        .seed = this->seed
                        }));

            history = std::vector<LatticeTrajectoryHistoryElement> ();
            history.reserve(this->history_chunk_size);
        }


        // increment step
        this->step++;

        // update_state
        lattice_network.update_state(state.lattice, std::ref(props), std::ref(this->state.homogeneous), next_reaction, 
                    event.site_one, event.site_two, latSolver.propensity_sum, latSolver.number_of_active_indices);


        // update_propensities 
        lattice_network.update_propensities(state.lattice, std::ref(this->state.homogeneous), this->update_function, 
                                        lattice_update_function, next_reaction, 
                                        event.site_one, event.site_two, props);


        return true;
    }
 
}

/* ---------------------------------------------------------------------------------------------- */
class NanoParticleSimulation : public Simulation<NanoSolver> {

    public:
    NanoParticle &nano_particle;
    std::vector<int> state;
    NanoSolver nanoSolver;
    std::vector<std::set<int>> site_reaction_dependency;

    std::vector<NanoTrajectoryHistoryElement> history;
    HistoryQueue<HistoryPacket<NanoTrajectoryHistoryElement>> &history_queue; 


    NanoParticleSimulation(NanoParticle &nano_particle,
               unsigned long int seed,
               int step,
               double time,
               std::vector<int> state,
               int history_chunk_size,
               HistoryQueue<HistoryPacket<NanoTrajectoryHistoryElement>> &history_queue
        ) :
        // call base class constructor
        Simulation<NanoSolver>(seed, history_chunk_size, step, time),
        nano_particle (nano_particle),
        state (state),
        history_queue(history_queue)
        { 
            history.reserve(this->history_chunk_size);
        };

    void init();
    bool execute_step();

};

void NanoParticleSimulation::init() {
    std::vector<std::set<int>> seed_site_reaction_dependency;
    std::vector<NanoReaction> seed_reactions;
    
    seed_site_reaction_dependency.resize(nano_particle.sites.size());
    nano_particle.compute_reactions(state, std::ref(seed_reactions), std::ref(seed_site_reaction_dependency));
    nanoSolver = NanoSolver(this->seed, std::ref(seed_reactions));
    site_reaction_dependency = seed_site_reaction_dependency;
    
}

bool NanoParticleSimulation::execute_step() {
    std::optional<Event> maybe_event = nanoSolver.event();

    if (! maybe_event) {

        return false;

    } else {
        // an event happens
        Event event = maybe_event.value();
        int next_reaction_id = event.index;
        NanoReaction next_reaction = nanoSolver.current_reactions[next_reaction_id];

        // update time
        this->time += event.dt;

        // record what happened
        history.push_back(NanoTrajectoryHistoryElement {
            .seed = this->seed,
            .reaction = next_reaction,
            .time = this->time,
            .step = this->step
            });

        if (history.size() == this->history_chunk_size ) {
            history_queue.insert_history(
                std::move(
                    HistoryPacket<NanoTrajectoryHistoryElement> {
                        .history = std::move(this->history),
                        .seed = this->seed
                        }));

            history = std::vector<NanoTrajectoryHistoryElement> ();
            history.reserve(this->history_chunk_size);
        }


        // increment step
        this->step++;

        // update state
        nano_particle.update_state(std::ref(state), next_reaction);

        // update list of current available reactions
        nano_particle.update_reactions(std::cref(state), next_reaction, std::ref(site_reaction_dependency), std::ref(nanoSolver.current_reactions));
        nanoSolver.update();

        return true;
    }
};


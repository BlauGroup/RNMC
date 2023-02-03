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

/*---------------------------------------------------------------------------*/

template <typename Solver, typename History>
class Simulation {
    public:
    unsigned long int seed;
    double time;
    int step; // number of reactions which have occoured
    unsigned long int history_chunk_size;
    std::vector<History> history;
    HistoryQueue<HistoryPacket<History>> &history_queue;   
    std::function<void(Update)> update_function;


    Simulation(unsigned long int seed,
               int history_chunk_size,
               HistoryQueue<HistoryPacket<History>> &history_queue, int step,
               double time
        ) :
        seed (seed),
        history_queue(history_queue),
        time (time),
        step (step),
        history_chunk_size (history_chunk_size)
        {
            history.reserve(history_chunk_size);
        };

    void execute_steps(int step_cutoff);
    void execute_time(double time_cutoff);

};

template <typename Solver, typename History>
void Simulation<Solver, History>::execute_steps(int step_cutoff) {

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

template <typename Solver, typename History>
void Simulation<Solver, History>::execute_time(double time_cutoff) {
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

template <typename Solver, typename History>
class ReactionNetworkSimulation : public Simulation<Solver, History> {
    private: 
        Solver solver;
    public:
    ReactionNetwork &GMC;
    std::vector<int> state;

    ReactionNetworkSimulation(ReactionNetwork &GMC, 
            unsigned long int seed,
            int history_chunk_size,
            HistoryQueue<HistoryPacket> &history_queue
        ) : 
        // time = 0, step = 0.0
        Simulation<LatSolver, History>(seed, history_chunk_size, history_queue, 0, 0.0),
        GMC (GMC),
        state (GMC.initial_state),
        { 
        };

    void init();
    bool execute_step();

};
template <typename Solver, typename History>
void ReactionNetworkSimulation<Solver, History>::init() {
    solver = Solver(seed, std::ref(GMC.initial_propensities));

}

template <typename Solver, typename History>
bool ReactionNetworkSimulation<Solver, History>::execute_step() {
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
        GMC.update_state(std::ref(state), next_reaction);


        // update propensities
        GMC.update_propensities(
            update_function,
            std::ref(state),
            next_reaction);

        return true;
    }
};

/* ---------------------------------------------------------------------- */
template<typename History>
class LatticeSimulation : public Simulation<LatSolver, History> {
    public:
    std::unordered_map<std::string,                     
        std::vector< std::pair<double, int> > > props;
    LatSolver latsolver;
    Lattice *lattice;
    std::function<void(LatticeUpdate, std::unordered_map<std::string,                     
                std::vector< std::pair<double, int> > > &)> lattice_update_function;
    std::function<void(Update)> update_function;

    LatticeSimulation(LGMC &LatticeReactionNetwork, unsigned long int seed, int history_chunk_size,
               HistoryQueue<HistoryPacket<History>> &history_queue) :
               // Call base class constructor, step = 0, time = 0.0
               Simulation<LatSolver, History>(seed, history_chunk_size, history_queue, 0, 0.0),
               latsolver (seed, std::ref(LGMC.initial_propensities)),
               lattice_update_function ([&] (LatticeUpdate lattice_update, 
               std::unordered_map<std::string,                     
                std::vector< std::pair<double, int> > > &props) {latsolver.update(lattice_update, props);}), 
               update_function ([&] (Update update) {latsolver.update(update);}) {
                
               };

    bool execute_step();
    void init(); 
    ~LatticeSimulation() { delete lattice;};

};

template<typename History>
void LatticeSimulation< History>::init() {
    
    lattice = new  Lattice(*LGMC.initial_lattice);
   
    LGMC.update_adsorp_state(this->lattice, this->props, latsolver.propensity_sum, latsolver.number_of_active_indices);
    LGMC.update_adsorp_props(this->lattice, lattice_update_function, this->state, std::ref(props));
    
}

template<typename History>
bool LatticeSimulation<History>::execute_step() {

    std::optional<LatticeEvent> maybe_event = latsolver.event_lattice(props);

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

        //if(site_1 == 5758) {
        //    std::cout << "next reaction: " << next_reaction << ", step: " << this->step << ", site 1: " << site_1 << ", site_2: " << site_2 << std::endl;
        //}
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
        this->LGMC.update_state(lattice, std::ref(props), std::ref(this->state), next_reaction, 
                    event.site_one, event.site_two, latsolver.propensity_sum, latsolver.number_of_active_indices);


        // update_propensities 
        this->LGMC.update_propensities(lattice, std::ref(this->state), this->update_function, 
                                        lattice_update_function, next_reaction, 
                                        event.site_one, event.site_two, props);


        return true;
    }
 
}

/* ---------------------------------------------------------------------------------------------- */

template <typename Solver, typename History>
class NanoParticleSimulation : public Simulation<Solver, History> {

    public:
    NanoParticle &NPMC;
    std::vector<int> state;
    Solver nanoSolver;
    HistoryQueue<StateHistoryPacket> &state_history_queue;
    std::vector<std::set<int>> site_reaction_dependency;

    NanoParticleSimulation(NanoParticle &NPMC,
               unsigned long int seed,
               int step,
               double time,
               std::vector<int> state,
               int history_chunk_size,
               HistoryQueue<HistoryPacket> &history_queue,
               HistoryQueue<StateHistoryPacket> &state_history_queue
        ) :
        // call base class constructor
        Simulation<LatSolver, History>(seed, history_chunk_size, history_queue, step, time),
        NPMC (NPMC),
        state (state),
        time (time),
        step (step),
        history_queue (history_queue),
        state_history_queue (state_history_queue),
        
        {
        };

    void init();
    bool execute_step();

};
template <typename History>
void NanoParticleSimulation<History>::init() {
    std::vector<std::set<int>> seed_site_reaction_dependency;
    std::vector<Reaction> seed_reactions;
    
    seed_site_reaction_dependency.resize(NPMC.sites.size());
    NPMC.compute_reactions(state, std::ref(seed_reactions), std::ref(seed_site_reaction_dependency));
    nanoSolver = NanoSolver(seed, std::ref(seed_reactions));
    site_reaction_dependency = seed_site_reaction_dependency;
    
}


template <typename History>
bool NanoParticleSimulation<History>::execute_step() {
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

template <typename History>
void NanoParticleSimulation::write_error_message(std::string s){
    char char_array[s.length()+1];
    strcpy(char_array, s.c_str());

    write(STDERR_FILENO, char_array, sizeof(char_array) - 1);
}

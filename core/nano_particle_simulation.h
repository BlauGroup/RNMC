#ifndef RNMC_NANO_PARTICLE_SIMULATION_H
#define RNMC_NANO_PARTICLE_SIMULATION_H

#include <vector>

#include "simulation.h"
#include "RNMC_types.h"
#include "queues.h"
#include "../NPMC/nano_solver.h"
#include "../NPMC/nano_particle.h"

class NanoParticleSimulation : public Simulation<NanoSolver> {
public:
    NanoParticle &nano_particle;
    std::vector<int> state;
    NanoSolver nanoSolver;
    std::vector<std::set<int>> site_reaction_dependency;
    std::vector<NanoTrajectoryHistoryElement> history;
    HistoryQueue<HistoryPacket<NanoTrajectoryHistoryElement>> &history_queue; 

    NanoParticleSimulation(NanoParticle &nano_particle,
               unsigned long int seed, int step,
               double time, std::vector<int> state,
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
    void print_output() {assert(true);};
};

#include "nano_particle_simulation.cpp"

#endif
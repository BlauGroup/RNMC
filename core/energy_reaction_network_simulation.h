#ifndef RNMC_ENERGY_REACTION_NETWORK_SIMULATION_H
#define RNMC_ENERGY_REACTION_NETWORK_SIMULATION_H

#include "../GMC/energy_reaction_network.h"
#include "simulation.h"

template <typename Solver>
class EnergyReactionNetworkSimulation : public Simulation<Solver> {
private: 
    Solver solver;
public:
    EnergyReactionNetwork &energy_reaction_network;
    EnergyState state;
    std::vector<ReactionNetworkTrajectoryHistoryElement> history;
    HistoryQueue<HistoryPacket<ReactionNetworkTrajectoryHistoryElement>> &history_queue; 

    EnergyReactionNetworkSimulation(EnergyReactionNetwork &reaction_network, 
            unsigned long int seed,
            int step,
            double time,
            EnergyState state,
            int history_chunk_size,
            HistoryQueue<HistoryPacket<ReactionNetworkTrajectoryHistoryElement>> &history_queue
        ) : 
        // call base class constructor
        Simulation<Solver>(seed, history_chunk_size, step, time),
        energy_reaction_network (reaction_network),
        state (state),
        history_queue(history_queue)
        { 
            history.reserve(this->history_chunk_size);
        };

    void init();
    bool execute_step();
    void print_output() {assert(true);};
};

#include "energy_reaction_network_simulation.cpp"

#endif 
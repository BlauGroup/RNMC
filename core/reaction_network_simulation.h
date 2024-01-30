#ifndef RNMC_REACTION_NETWORK_SIMULATION_H
#define RNMC_REACTION_NETWORK_SIMULATION_H

#include "../GMC/gillespie_reaction_network.h"
#include "simulation.h"

template <typename Solver>
class ReactionNetworkSimulation : public Simulation<Solver> {
private: 
    Solver solver;
public:
    GillespieReactionNetwork &reaction_network;
    std::vector<int> state;
    std::vector<ReactionNetworkTrajectoryHistoryElement> history;
    HistoryQueue<HistoryPacket<ReactionNetworkTrajectoryHistoryElement>> &history_queue; 

    ReactionNetworkSimulation(GillespieReactionNetwork reaction_network, 
            unsigned long int seed,
            int step,
            double time,
            std::vector<int> state,
            int history_chunk_size,
            HistoryQueue<HistoryPacket<ReactionNetworkTrajectoryHistoryElement>> &history_queue
        ) : 
        // call base class constructor
        Simulation<Solver>(seed, history_chunk_size, step, time),
        reaction_network (reaction_network),
        state (state),
        history_queue(history_queue)
        { 
            history.reserve(this->history_chunk_size);
        };

    void init();
    bool execute_step();
    void print_output() {assert(true);};
};

#include "reaction_network_simulation.cpp"

#endif 
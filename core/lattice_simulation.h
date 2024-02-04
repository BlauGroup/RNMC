#ifndef RNMC_LATTICE_SIMULATION_H
#define RNMC_LATTICE_SIMULATION_H

#include <map>
#include <vector>
#include <memory>

#include "dispatcher.h"
#include "simulation.h"
#include "../LGMC/lattice_solver.h"
#include "../LGMC/lattice_reaction_network.h"

class LatticeSimulation : public Simulation<LatticeSolver> {
    public:
    std::unordered_map<std::string, std::vector< std::pair<double, int> > > props;
    LatticeSolver latSolver;
    LatticeReactionNetwork lattice_network;
    LatticeState state;
    std::function<void(LatticeUpdate, std::unordered_map<std::string,                     
                        std::vector< std::pair<double, int> > > &)> 
                        lattice_update_function;
    
    std::function<void(Update)> update_function;
    std::vector<LatticeTrajectoryHistoryElement> history;
    HistoryQueue<HistoryPacket<LatticeTrajectoryHistoryElement>> &history_queue; 

    LatticeSimulation(LatticeReactionNetwork lattice_network, unsigned long int seed, 
                      int step, double time, LatticeState state_in, int history_chunk_size,
                      HistoryQueue<HistoryPacket<LatticeTrajectoryHistoryElement>> &history_queue
        ) :
        // call base class constructor
        Simulation<LatticeSolver>(seed, history_chunk_size, step, time),
        lattice_network(lattice_network),
        history_queue(history_queue)
        { 
            state.homogeneous = state_in.homogeneous;
            state.lattice = std::move(state_in.lattice);
            history.reserve(this->history_chunk_size);
        };
    
    void init(); 
    bool execute_step();
    void print_output();
    // ~LatticeSimulation() { delete state.lattice;};
};

#include "lattice_simulation.cpp"
#endif
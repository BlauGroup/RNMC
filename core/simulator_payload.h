#ifndef RNMC_SIMULATOR_PAYLOAD_H
#define RNMC_SIMULATOR_PAYLOAD_H

#include "simulation.h"
#include "RNMC_types.h"
#include "queues.h"

/* ---------------------------------------------------------------------- 
    size of history chunks which we write to the database.
    if you make this too small, it will force the dispatcher to
    perform lots of really small DB transactions which is bad.
    20000 is a good value. Only change this if you fully understand the
    performance implications 
 ---------------------------------------------------------------------- */

constexpr int history_chunk_size = 20000;

template <
    typename Solver, 
    typename Model, 
    typename StateHistory, 
    typename TrajHistory, 
    typename CutoffHistory, 
    typename Sim, 
    typename State>

class SimulatorPayload {

public: 
    Model &model;
    HistoryQueue<HistoryPacket<TrajHistory>> &history_queue;
    HistoryQueue<HistoryPacket<StateHistory>> &state_history_queue;
    HistoryQueue<HistoryPacket<CutoffHistory>> &cutoff_history_queue;
    SeedQueue &seed_queue;
    Cutoff cutoff;
    std::vector<bool>::iterator running;
    std::map<int, State> seed_state_map;
    std::map<int, int> seed_step_map;
    std::map<int, double> seed_time_map;

    SimulatorPayload(
        Model &model,
        HistoryQueue<HistoryPacket<TrajHistory>> &history_queue,
        HistoryQueue<HistoryPacket<StateHistory>> &state_history_queue,
        HistoryQueue<HistoryPacket<CutoffHistory>> &cutoff_history_queue,
        SeedQueue &seed_queue,
        Cutoff cutoff,
        std::vector<bool>::iterator running, 
        std::map<int, State> seed_state_map,
        std::map<int, int> seed_step_map,
        std::map<int, double> seed_time_map
        ):
            model (model),
            history_queue (history_queue),
            state_history_queue (state_history_queue),
            cutoff_history_queue (cutoff_history_queue),
            seed_queue (seed_queue),
            cutoff (cutoff),
            running (running),
            seed_state_map (seed_state_map),
            seed_step_map (seed_step_map),
            seed_time_map (seed_time_map)
        {};

    void run_simulator() {

          while (std::optional<unsigned long int> maybe_seed =
               seed_queue.get_seed()) {
            
            unsigned long int seed = maybe_seed.value();
            int step = seed_step_map[seed];
            double time = seed_time_map[seed];
            State state = seed_state_map[seed];

            Sim simulation(model, seed, step, time, state, 
                           history_chunk_size, history_queue);
            simulation.init();

            switch(cutoff.type_of_cutoff) {
            case step_termination :
                simulation.execute_steps(cutoff.bound.step);
                break;
            case time_termination :
                simulation.execute_time(cutoff.bound.time);
                break;
            }

            simulation.print_output();
            
            // Make a vector of StateHistoryElements for the current state
            std::vector<StateHistory> state_packet;

            // Make a vector of CutoffHistories for the current time and step
            std::vector<CutoffHistory> cutoff_packet;

            model.store_checkpoint(state_packet, simulation.state, seed, 
                                 simulation.step, simulation.time, cutoff_packet);

            // Construct a history packet from the history elements and add it to the queue
            state_history_queue.insert_history(
                std::move(
                    HistoryPacket<StateHistory> {
                        .seed = seed,
                        .history = state_packet
                    }));

            // Construct a history packet from the history elements and add it to the queue
            cutoff_history_queue.insert_history(
                std::move(
                    HistoryPacket<CutoffHistory> {
                        .seed = seed,
                        .history = cutoff_packet
                    }));

            // Move the remainder of the history into the queue to be saved
            history_queue.insert_history(
                std::move(
                    HistoryPacket<TrajHistory> {
                        .history = std::move(simulation.history),
                        .seed = seed
                        }));

        }

        *running = false;
    };
};

#endif 
#pragma once
#include <mutex>
#include <thread>
#include "simulation.h"
#include "../core/queues.h"

struct GMCHistoryPacket {
    std::vector<HistoryElement> history;
    unsigned long int seed;
};



template <typename Solver, typename Model>
struct SimulatorPayload {
    ReactionNetwork &reaction_network;
    HistoryQueue<GMCHistoryPacket> &history_queue;
    SeedQueue &seed_queue;
    int step_cutoff;

    SimulatorPayload(
        ReactionNetwork &reaction_network,
        HistoryQueue<GMCHistoryPacket> &history_queue,
        SeedQueue &seed_queue,
        int step_cutoff
        ) :

            reaction_network (reaction_network),
            history_queue (history_queue),
            seed_queue (seed_queue),
            step_cutoff (step_cutoff) {};

    void run_simulator() {

        while (std::optional<unsigned long int> maybe_seed =
               seed_queue.get_seed()) {

            unsigned long int seed = maybe_seed.value();
            Simulation<Solver, Model> simulation (reaction_network, seed, step_cutoff);
            simulation.execute_steps(step_cutoff);

            // Calling resize() with a smaller size has no effect on the capacity of a vector.
            // It will not free memory.
            simulation.history.resize(simulation.step);
            history_queue.insert_history(
                std::move(
                    GMCHistoryPacket {
                        .history = std::move(simulation.history),
                        .seed = seed
                        }));
        }
    }
};

template <typename Solver, typename Model>
struct Dispatcher {
    SqlConnection reaction_database;
    SqlConnection initial_state_database;
    ReactionNetwork reaction_network;
    SqlStatement<TrajectoriesSql> trajectories_stmt;
    SqlWriter<TrajectoriesSql> trajectories_writer;
    HistoryQueue<GMCHistoryPacket> history_queue;
    SeedQueue seed_queue;
    std::vector<std::thread> threads;
    int step_cutoff;
    int number_of_simulations;
    int number_of_threads;

    Dispatcher(
        std::string reaction_database_file,
        std::string initial_state_database_file,
        unsigned long int number_of_simulations,
        unsigned long int base_seed,
        int number_of_threads,
        int step_cutoff,
        int dependency_threshold) :
        reaction_database (
            reaction_database_file,
            SQLITE_OPEN_READWRITE),
        initial_state_database (
            initial_state_database_file,
            SQLITE_OPEN_READWRITE),
        reaction_network (
            reaction_database,
            initial_state_database,
            dependency_threshold),
        trajectories_stmt (initial_state_database),
        trajectories_writer (trajectories_stmt),
        history_queue (),
        seed_queue (number_of_simulations, base_seed),

        // don't want to start threads in the constructor.
        threads (),
        step_cutoff (step_cutoff),
        number_of_simulations (number_of_simulations),
        number_of_threads (number_of_threads)
        {};

    void run_dispatcher();
    void record_simulation_history(GMCHistoryPacket history_packet);
};



template <typename Solver, typename Model>
void Dispatcher<Solver, Model>::run_dispatcher() {

    threads.resize(number_of_threads);
    for (int i = 0; i < number_of_threads; i++) {
        threads[i] = std::thread (
            [](SimulatorPayload<Solver, Model> payload) {payload.run_simulator();},
            SimulatorPayload<Solver, Model> (
                reaction_network,
                history_queue,
                seed_queue,
                step_cutoff)
            );

    }

    int trajectories_written = 0;
    while (trajectories_written < number_of_simulations) {

        std::optional<GMCHistoryPacket>
            maybe_history_packet = history_queue.get_history();

        if (maybe_history_packet) {
            GMCHistoryPacket history_packet = std::move(maybe_history_packet.value());
            record_simulation_history(std::move(history_packet));
            trajectories_written += 1;
        };
    }

    for (int i = 0; i < number_of_threads; i++) threads[i].join();

    initial_state_database.exec(
        "DELETE FROM trajectories WHERE rowid NOT IN"
        "(SELECT MIN(rowid) FROM trajectories GROUP BY seed, step);");

    std::cerr << time_stamp()
              << "removing duplicate trajectories...\n";


};

template <typename Solver, typename Model>
void Dispatcher<Solver, Model>::record_simulation_history(GMCHistoryPacket history_packet) {
    int count = 0;
    constexpr int transaction_size = 20000;
    initial_state_database.exec("BEGIN");
    for (unsigned long int i = 0; i < history_packet.history.size(); i++) {
        trajectories_writer.insert(

            TrajectoriesSql {
                .seed = (int) history_packet.seed,
                .step = (int) i,
                .reaction_id = history_packet.history[i].reaction_id,
                .time = history_packet.history[i].time
            });
        count++;
        if (count % transaction_size == 0) {
            initial_state_database.exec("COMMIT;");
            initial_state_database.exec("BEGIN");

        }
    }
    initial_state_database.exec("COMMIT;");

    std::cerr << time_stamp()
              << "wrote trajectory "
              << history_packet.seed
              << " to database\n";
};

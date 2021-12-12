#pragma once
#include <mutex>
#include <thread>
#include "sql.h"
#include "simulation.h"
#include "queues.h"

struct HistoryPacket {
    std::vector<HistoryElement> history;
    unsigned long int seed;
};



template <typename Solver, typename Model>
struct SimulatorPayload {
    Model &model;
    HistoryQueue<HistoryPacket> &history_queue;
    SeedQueue &seed_queue;
    int step_cutoff;

    SimulatorPayload(
        Model &model,
        HistoryQueue<HistoryPacket> &history_queue,
        SeedQueue &seed_queue,
        int step_cutoff
        ) :

            model (model),
            history_queue (history_queue),
            seed_queue (seed_queue),
            step_cutoff (step_cutoff) {};

    void run_simulator() {

        while (std::optional<unsigned long int> maybe_seed =
               seed_queue.get_seed()) {

            unsigned long int seed = maybe_seed.value();
            Simulation<Solver, Model> simulation (model, seed, step_cutoff);
            simulation.execute_steps(step_cutoff);

            // Calling resize() with a smaller size has no effect on the capacity of a vector.
            // It will not free memory.
            simulation.history.resize(simulation.step);
            history_queue.insert_history(
                std::move(
                    HistoryPacket {
                        .history = std::move(simulation.history),
                        .seed = seed
                        }));
        }
    }
};

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename TrajectoriesSql>

struct Dispatcher {
    SqlConnection model_database;
    SqlConnection initial_state_database;
    Model model;
    SqlStatement<TrajectoriesSql> trajectories_stmt;
    SqlWriter<TrajectoriesSql> trajectories_writer;
    HistoryQueue<HistoryPacket> history_queue;
    SeedQueue seed_queue;
    std::vector<std::thread> threads;
    int step_cutoff;
    int number_of_simulations;
    int number_of_threads;

    Dispatcher(
        std::string model_database_file,
        std::string initial_state_database_file,
        unsigned long int number_of_simulations,
        unsigned long int base_seed,
        int number_of_threads,
        int step_cutoff,
        Parameters parameters) :
        model_database (
            model_database_file,
            SQLITE_OPEN_READWRITE),
        initial_state_database (
            initial_state_database_file,
            SQLITE_OPEN_READWRITE),
        model (
            model_database,
            initial_state_database,
            parameters),
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
    void record_simulation_history(HistoryPacket history_packet);
};



template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename TrajectoriesSql>

void Dispatcher<Solver, Model, Parameters, TrajectoriesSql>::run_dispatcher() {

    threads.resize(number_of_threads);
    for (int i = 0; i < number_of_threads; i++) {
        threads[i] = std::thread (
            [](SimulatorPayload<Solver, Model> payload) {payload.run_simulator();},
            SimulatorPayload<Solver, Model> (
                model,
                history_queue,
                seed_queue,
                step_cutoff)
            );

    }

    int trajectories_written = 0;
    while (trajectories_written < number_of_simulations) {

        std::optional<HistoryPacket>
            maybe_history_packet = history_queue.get_history();

        if (maybe_history_packet) {
            HistoryPacket history_packet = std::move(maybe_history_packet.value());
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

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename TrajectoriesSql
    >
void Dispatcher<Solver, Model, Parameters, TrajectoriesSql>::record_simulation_history(HistoryPacket history_packet) {
    int count = 0;
    constexpr int transaction_size = 20000;
    initial_state_database.exec("BEGIN");
    for (unsigned long int i = 0; i < history_packet.history.size(); i++) {
        trajectories_writer.insert(
            model.history_element_to_sql(
                (int) history_packet.seed,
                (int) i,
                history_packet.history[i]));
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

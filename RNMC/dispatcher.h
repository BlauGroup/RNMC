#pragma once
#include <queue>
#include <mutex>
#include <thread>
#include "simulation.h"


struct SeedQueue {
    std::queue<unsigned long int> seeds;
    std::mutex mutex;

    SeedQueue(unsigned long int number_of_seeds, unsigned long int base_seed) {
        for (unsigned long int i = base_seed;
             i < number_of_seeds + base_seed;
             i++) {
            seeds.push(i);
        }
    }

    std::optional<unsigned long int> get_seed() {
        std::lock_guard<std::mutex> lock (mutex);

        if (seeds.empty()) {
            return std::optional<unsigned long int> ();
        } else {
            unsigned long int result = seeds.front();
            seeds.pop();
            return std::optional<unsigned long int> (result);
        }
    }
};


struct HistoryQueue {
    std::queue<std::vector<HistoryElement>> histories;
    std::mutex mutex;

    void insert_history(std::vector<HistoryElement> &&history) {
        std::lock_guard<std::mutex> lock (mutex);
        histories.push(std::move(history));
    }

    std::vector<HistoryElement> get_history() {
        std::lock_guard<std::mutex> lock (mutex);
        std::vector<HistoryElement> result = std::move(histories.front());
        histories.pop();
        return result;
    };

};


template <typename Solver>
struct SimulatorPayload {
    ReactionNetwork &reaction_network;
    HistoryQueue &history_queue;
    SeedQueue &seed_queue;
    int step_cutoff;
    bool &running;

    SimulatorPayload(
        ReactionNetwork &reaction_network,
        HistoryQueue &history_queue,
        SeedQueue &seed_queue,
        int step_cutoff,
        bool &running
        ) :

            reaction_network (reaction_network),
            history_queue (history_queue),
            seed_queue (seed_queue),
            step_cutoff (step_cutoff),
            running (running) {};

    void run_simulator() {

        while (std::optional<unsigned long int> maybe_seed =
               seed_queue.get_seed()) {

            unsigned long int seed = maybe_seed.value();
            Simulation<Solver> simulation (reaction_network, seed, step_cutoff);
            simulation.execute_steps(step_cutoff);
            history_queue.insert_history(std::move(simulation.history));
        }

        running = false;

    }
};


struct Dispatcher {
    SqlConnection reaction_database;
    SqlConnection initial_state_database;
    ReactionNetwork reaction_network;
    SqlStatement<TrajectoriesSql> trajectories_stmt;
    SqlWriter<TrajectoriesSql> trajectories_writer;
    HistoryQueue history_queue;
    SeedQueue seed_queue;
    std::vector<std::thread> threads;
    std::vector<bool> running;
    int step_cutoff;

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
        running (number_of_threads, true),
        step_cutoff (step_cutoff)
        {};
};


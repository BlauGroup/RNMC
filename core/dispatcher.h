#pragma once
#include <mutex>
#include <thread>
#include "sql.h"
#include "simulation.h"
#include "queues.h"
#include <csignal>
#include <iostream>
#include <string>



enum TypeOfCutoff { step_termination, time_termination };

struct Cutoff {
    union  { int step; double time; } bound;
    TypeOfCutoff type_of_cutoff;
};


// size of history chunks which we write to the DB.
// if you make this too small, it will force the dispatcher to
// perform lots of really small DB transactions which is bad.
// 20000 is a good value. Only change this if you fully understand the
// performance implications

constexpr int history_chunk_size = 100;

template <typename Solver, typename Model, typename History>
struct SimulatorPayload {
    Model &model;
    HistoryQueue<HistoryPacket<History>> &history_queue;
    SeedQueue &seed_queue;
    Cutoff cutoff;
    std::vector<bool>::iterator running;
    bool isLGMC;

    SimulatorPayload(
        Model &model,
        HistoryQueue<HistoryPacket<History>> &history_queue,
        SeedQueue &seed_queue,
        Cutoff cutoff,
        std::vector<bool>::iterator running, 
        bool isLGMC
        ) :

            model (model),
            history_queue (history_queue),
            seed_queue (seed_queue),
            cutoff (cutoff),
            running (running),
            isLGMC (isLGMC)
        {};

    void run_simulator() {

        while (std::optional<unsigned long int> maybe_seed =
               seed_queue.get_seed()) {

            unsigned long int seed = maybe_seed.value();

            if(!isLGMC) {
                Simulation<Solver, Model, History> simulation (model, seed, history_chunk_size, history_queue);
                simulation.init(model);

                switch(cutoff.type_of_cutoff) {
                case step_termination :
                    simulation.execute_steps(cutoff.bound.step);
                    break;
                case time_termination :
                    simulation.execute_time(cutoff.bound.time);
                    break;
                }
             
                    history_queue.insert_history(
                    std::move(
                    HistoryPacket<History> {
                    .history = std::move(simulation.history),
                    .seed = seed
                    }));

            }
            else {
                LatticeSimulation<Model, History> simulation (model, seed, history_chunk_size, history_queue);
                simulation.init();

                switch(cutoff.type_of_cutoff) {
                case step_termination :
                    simulation.execute_steps(cutoff.bound.step);
                    break;
                case time_termination :
                    simulation.execute_time(cutoff.bound.time);
                    break;
                }

                    history_queue.insert_history(
                    std::move(
                    HistoryPacket<History> {
                    .history = std::move(simulation.history),
                    .seed = seed
                    }));
            }


        }

        *running = false;
    }
};



template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename OutputSQL, 
    typename History>

struct Dispatcher {
    SqlConnection model_database;
    SqlConnection initial_state_database;
    Model model;
    SqlStatement<OutputSQL> trajectories_stmt;
    SqlWriter<OutputSQL> trajectories_writer;
    HistoryQueue<HistoryPacket<History>> history_queue;
    SeedQueue seed_queue;
    std::vector<std::thread> threads;
    std::vector<bool> running;
    Cutoff cutoff;
    TypeOfCutoff type_of_cutoff;
    int number_of_simulations;
    int number_of_threads;
    bool isLGMC;

    
    Dispatcher(
        std::string model_database_file,
        std::string initial_state_database_file,
        unsigned long int number_of_simulations,
        unsigned long int base_seed,
        int number_of_threads,
        Cutoff cutoff,
        Parameters parameters, 
        bool isLGMC) :
        model_database (
            model_database_file,
            SQLITE_OPEN_READWRITE),
        initial_state_database (
            initial_state_database_file,
            SQLITE_OPEN_READWRITE),
        model (model_database,
            initial_state_database,
            parameters),
        trajectories_stmt (initial_state_database),
        trajectories_writer (trajectories_stmt),
        history_queue (),
        seed_queue (number_of_simulations, base_seed),

        // don't want to start threads in the constructor.
        threads (),
        running (number_of_threads, false),
        cutoff (cutoff),
        number_of_simulations (number_of_simulations),
        number_of_threads (number_of_threads), 
        isLGMC(isLGMC)
        {
        };

    void run_dispatcher();
    void record_simulation_history(HistoryPacket<History> history_packet);
};

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename OutputSQL, 
    typename History>

void Dispatcher<Solver, Model, Parameters, OutputSQL, History>::run_dispatcher() {


    threads.resize(number_of_threads);
    for (int i = 0; i < number_of_threads; i++) {
        running[i] = true;
        threads[i] = std::thread (
            [](SimulatorPayload<Solver, Model, History> payload) {payload.run_simulator();},
            SimulatorPayload<Solver, Model, History> (
                model,
                history_queue,
                seed_queue,
                cutoff,
                running.begin() + i, 
                isLGMC
                )
            );

    }

    bool finished = false;
    while (! finished) {

        if ( history_queue.empty() ) {

            bool all_simulations_finished = true;

            for ( bool flag : running ) {
                if ( flag ) {
                    all_simulations_finished = false;
                    break;
                }
            }

            if ( all_simulations_finished )
                finished = true;
        }

        std::optional<HistoryPacket<History>>
        maybe_history_packet = history_queue.get_history();

        if (maybe_history_packet) {
            HistoryPacket<History> history_packet = std::move(maybe_history_packet.value());
            record_simulation_history(std::move(history_packet));
        }

        
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
    typename OutputSQL,
    typename History>
void Dispatcher<Solver, Model, Parameters, OutputSQL, History>::record_simulation_history(HistoryPacket<History> history_packet) {
    initial_state_database.exec("BEGIN");


    for (unsigned long int i = 0; i < history_packet.history.size(); i++) {
        trajectories_writer.insert(
            model.history_element_to_sql(
                (int) history_packet.seed,
                history_packet.history[i]));
    }
    initial_state_database.exec("COMMIT;");

    std::cerr << time_stamp()
              << "wrote "
              << history_packet.history.size()
              << " events from trajectory "
              << history_packet.seed
              << " to database\n";
};

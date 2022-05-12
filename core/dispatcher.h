#pragma once
#include <mutex>
#include <thread>
#include "sql.h"
#include "simulation.h"
#include "queues.h"
#include <csignal>



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

constexpr int history_chunk_size = 20000;

template <typename Solver, typename Model>
struct SimulatorPayload {
    Model &model;
    HistoryQueue<HistoryPacket> &history_queue;
    HistoryQueue<StateHistoryPacket> &state_history_queue;
    SeedQueue &seed_queue;
    Cutoff cutoff;
    std::vector<bool>::iterator running;

    SimulatorPayload(
        Model &model,
        HistoryQueue<HistoryPacket> &history_queue,
        HistoryQueue<StateHistoryPacket> &state_history_queue,
        SeedQueue &seed_queue,
        Cutoff cutoff,
        std::vector<bool>::iterator running
        ) :

            model (model),
            history_queue (history_queue),
            state_history_queue (state_history_queue),
            seed_queue (seed_queue),
            cutoff (cutoff),
            running (running)
        {};

    void run_simulator() {

        while (std::optional<unsigned long int> maybe_seed =
               seed_queue.get_seed()) {

            unsigned long int seed = maybe_seed.value();

            // Generate the model

            Simulation<Solver, Model> simulation (model, seed, history_chunk_size, history_queue, state_history_queue);


            switch(cutoff.type_of_cutoff) {
            case step_termination :
                simulation.execute_steps(cutoff.bound.step);
                break;
            case time_termination :
                simulation.execute_time(cutoff.bound.time);
                break;
            }

            // Make a vector of StateHistoryElements for the current state
            std::vector<StateHistoryElement> state_packet;
            for (unsigned int i = 0; i < model.sites.size(); i++) {
                state_packet.push_back(StateHistoryElement{
                    .seed = seed,
                    .site_id = i,
                    .degree_of_freedom = simulation.state[i]
                });
            }
            // Construct a history packet from the history elements and add it to the queue
            state_history_queue.insert_history(
                std::move(
                    StateHistoryPacket {
                        .seed = seed,
                        .state = state_packet
                    }));

            // Move the remainder of the history into the queue to be saved
            history_queue.insert_history(
                std::move(
                    HistoryPacket {
                        .history = std::move(simulation.history),
                        .seed = seed
                        }));

        }

        *running = false;
    }
};



template <
    typename Solver,
    typename Model,
    typename WriteTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql>

struct Dispatcher {
    SqlConnection model_database;
    SqlConnection initial_state_database;
    Model model;
    SqlStatement<WriteTrajectoriesSql> trajectories_stmt;
    SqlWriter<WriteTrajectoriesSql> trajectories_writer;
    SqlStatement<WriteStateSql> state_stmt;
    SqlWriter<WriteStateSql> state_writer;
    HistoryQueue<HistoryPacket> history_queue;
    HistoryQueue<StateHistoryPacket> state_history_queue;
    SeedQueue seed_queue;
    std::vector<std::thread> threads;
    std::vector<bool> running;
    Cutoff cutoff;
    TypeOfCutoff type_of_cutoff;
    int number_of_simulations;
    int number_of_threads;
    struct sigaction action;

    Dispatcher(
        std::string model_database_file,
        std::string initial_state_database_file,
        unsigned long int number_of_simulations,
        unsigned long int base_seed,
        int number_of_threads,
        Cutoff cutoff) :
        model_database (
            model_database_file,
            SQLITE_OPEN_READWRITE),
        initial_state_database (
            initial_state_database_file,
            SQLITE_OPEN_READWRITE),
        model (
            model_database,
            initial_state_database),
        trajectories_stmt (initial_state_database),
        trajectories_writer (trajectories_stmt),
        state_stmt (initial_state_database),
        state_writer (state_stmt),
        history_queue (),
        state_history_queue (),
        seed_queue (number_of_simulations, base_seed),

        // don't want to start threads in the constructor.
        threads (),
        running (number_of_threads, false),
        cutoff (cutoff),
        number_of_simulations (number_of_simulations),
        number_of_threads (number_of_threads)
        {
        };

    void run_dispatcher();
    void record_simulation_history(HistoryPacket history_packet);
    void record_state(StateHistoryPacket state_history_packet);
};

void signalHandler(int signum) {
    
    write_error_message("received signal " + std::to_string(signum) + "\n");
    switch (signum) {
        case 2:
        case 4:
        case 6:
            // A SIGTERM has been issued, handle this gracefully 
            // by writing current states and trajectories to the db

            write_error_message("SIGTERM received. Terminating NPMC run(s) early.\n");
            do_shutdown = 1;
            shutdown_requested = true;

            write_error_message("Writing current states and history queue to the initial_states db\n");
            break;
        default: exit(signum);
    }
    
}

template <
    typename Solver,
    typename Model,
    typename WriteTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql>

void Dispatcher<Solver, Model, WriteTrajectoriesSql, WriteStateSql, ReadStateSql>::run_dispatcher() {

    // Check if there are existing runs
        // Check if individual states have been written to the initial_state.sqlite

        // Check if there are any trajectories


    // Create a set of masks containing the SIGTERM.
    sigset_t mask;
    sigemptyset (&mask);
    sigaddset (&mask, SIGCONT);
    sigaddset (&mask, SIGTERM);
    
    // Set the masks so the child threads inherit the sigmask to ignore SIGTERM
    pthread_sigmask(SIG_BLOCK, &mask, NULL);

    threads.resize(number_of_threads);
    for (int i = 0; i < number_of_threads; i++) {
        running[i] = true;
        threads[i] = std::thread (
            [](SimulatorPayload<Solver, Model> payload) {payload.run_simulator();},
            SimulatorPayload<Solver, Model> (
                model,
                history_queue,
                state_history_queue,
                seed_queue,
                cutoff,
                running.begin() + i
                )
            );

    }
    // Unset the sigmask so that the parent thread resumes catching errors as normal
    pthread_sigmask(SIG_UNBLOCK, &mask, NULL);
    
    action.sa_handler = signalHandler;
    sigemptyset(&action.sa_mask);
    action.sa_flags = 0;
    sigaction(SIGINT, &action, NULL);

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


        std::optional<HistoryPacket>
            maybe_history_packet = history_queue.get_history();

        if (maybe_history_packet) {
            HistoryPacket history_packet = std::move(maybe_history_packet.value());
            record_simulation_history(std::move(history_packet));
        };

        std::optional<StateHistoryPacket>
            maybe_state_history_packet = state_history_queue.get_history();

        if (maybe_state_history_packet) {
            StateHistoryPacket state_history_packet = std::move(maybe_state_history_packet.value());
            record_state(std::move(state_history_packet));
        };
    }

    for (int i = 0; i < number_of_threads; i++) threads[i].join();

    initial_state_database.exec(
        "DELETE FROM trajectories WHERE rowid NOT IN"
        "(SELECT MIN(rowid) FROM trajectories GROUP BY seed, step);");

    initial_state_database.close();
    model_database.close();
    std::cerr << time_stamp()
              << "removing duplicate trajectories...\n";


};

template <
    typename Solver,
    typename Model,
    typename WriteTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql
    >
void Dispatcher<Solver, Model, WriteTrajectoriesSql, WriteStateSql, ReadStateSql>::record_simulation_history(HistoryPacket history_packet) {
    initial_state_database.exec("BEGIN;");


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

template <
    typename Solver,
    typename Model,
    typename WriteTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql
    >
void Dispatcher<Solver, Model, WriteTrajectoriesSql, WriteStateSql, ReadStateSql>::record_state(StateHistoryPacket state_history_packet) {

    // Wipe the database of the states corresponding to this seed. This gets rid of the previously written states
    std::string delete_statement = "DELETE FROM interupt_state WHERE seed = " + std::to_string(state_history_packet.seed) + ";";
    
    initial_state_database.exec(delete_statement);
    
    initial_state_database.exec("BEGIN;");

    for (unsigned long int i = 0; i < state_history_packet.state.size(); i++) {
        state_writer.insert(
            model.state_history_element_to_sql(
                (int) state_history_packet.seed,
                state_history_packet.state[i]));
    }
    initial_state_database.exec("COMMIT;");

    std::cerr << time_stamp()
              << "wrote "
              << state_history_packet.state.size()
              << " states for trajectory "
              << state_history_packet.seed
              << " to database\n";
};

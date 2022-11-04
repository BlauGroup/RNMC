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
    HistoryQueue<CutoffHistoryPacket> &cutoff_history_queue;
    SeedQueue &seed_queue;
    Cutoff cutoff;
    std::vector<bool>::iterator running;
    std::map<int, std::vector<int>> seed_state_map;
    std::map<int, int> seed_step_map;
    std::map<int, double> seed_time_map;

    SimulatorPayload(
        Model &model,
        HistoryQueue<HistoryPacket> &history_queue,
        HistoryQueue<StateHistoryPacket> &state_history_queue,
        HistoryQueue<CutoffHistoryPacket> &cutoff_history_queue,
        SeedQueue &seed_queue,
        Cutoff cutoff,
        std::vector<bool>::iterator running,
        std::map<int, std::vector<int>> seed_state_map,
        std::map<int, int> seed_step_map,
        std::map<int, double> seed_time_map
        ) :

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
        // int seed = 1000;
        // ReadStateSql sql_reader;
        while (std::optional<unsigned long int> maybe_seed =
               seed_queue.get_seed()) {

            unsigned long int seed = maybe_seed.value();

            int step = seed_step_map[seed];
            double time = seed_time_map[seed];
            std::vector<int> state = seed_state_map[seed];
            std::vector<Reaction> seed_reactions;
            std::vector<std::set<int>> seed_site_reaction_dependency;
            seed_site_reaction_dependency.resize(model.sites.size());
            model.compute_reactions(state, std::ref(seed_reactions), std::ref(seed_site_reaction_dependency));

            Simulation<Solver, Model> simulation (model, seed, step, time, state, seed_reactions, seed_site_reaction_dependency, history_chunk_size, history_queue, state_history_queue);


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
                    .site_id = static_cast<int>(i),
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

            // Make a vector of CutoffHistoryElements for the current time and step
            std::vector<CutoffHistoryElement> cutoff_packet;

            cutoff_packet.push_back(CutoffHistoryElement {
                .seed = seed,
                .step = simulation.step,
                .time = simulation.time
            });

            // Construct a history packet from the history elements and add it to the queue
            cutoff_history_queue.insert_history(
                std::move(
                    CutoffHistoryPacket {
                        .seed = seed,
                        .cutoff = cutoff_packet
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
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql>

struct Dispatcher {
    SqlConnection model_database;
    SqlConnection initial_state_database;
    Model model;
    SqlStatement<WriteTrajectoriesSql> trajectories_stmt;
    SqlWriter<WriteTrajectoriesSql> trajectories_writer;
    SqlStatement<WriteStateSql> state_stmt;
    SqlWriter<WriteStateSql> state_writer;
    SqlStatement<WriteCutoffSql> cutoff_stmt;
    SqlWriter<WriteCutoffSql> cutoff_writer;
    HistoryQueue<HistoryPacket> history_queue;
    HistoryQueue<StateHistoryPacket> state_history_queue;
    HistoryQueue<CutoffHistoryPacket> cutoff_history_queue;
    SeedQueue seed_queue;
    std::vector<std::thread> threads;
    std::vector<bool> running;
    Cutoff cutoff;
    TypeOfCutoff type_of_cutoff;
    int number_of_simulations;
    int number_of_threads;
    struct sigaction action;
    std::map<int, std::vector<int>> seed_state_map;
    std::map<int, int> seed_step_map;
    std::map<int, double> seed_time_map;

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
        cutoff_stmt (initial_state_database),
        cutoff_writer (cutoff_stmt),
        history_queue (),
        state_history_queue (),
        cutoff_history_queue (),
        seed_queue (number_of_simulations, base_seed),

        // don't want to start threads in the constructor.
        threads (),
        running (number_of_threads, false),
        cutoff (cutoff),
        number_of_simulations (number_of_simulations),
        number_of_threads (number_of_threads),
        seed_state_map (),
        seed_step_map (),
        seed_time_map ()
        {
            std::vector<int> default_state = model.initial_state;

            SeedQueue temp_seed_queue = SeedQueue(number_of_simulations, base_seed);
            std::map<int, std::vector<int>> temp_seed_state_map;
            std::map<int, int> temp_seed_step_map;
            std::map<int, double> temp_seed_time_map;

            while (std::optional<unsigned long int> maybe_seed =
               temp_seed_queue.get_seed()){
                unsigned long int seed = maybe_seed.value();
                temp_seed_state_map.insert(std::make_pair(seed, default_state));
                temp_seed_step_map.insert(std::make_pair(seed, 0));
                temp_seed_time_map.insert(std::make_pair(seed, 0.0));
            }

            // Populate the seed-state map
            SqlStatement<ReadStateSql> state_statement(initial_state_database);
            SqlReader<ReadStateSql> state_reader(state_statement);
            bool read_interupt_states = false;
            while (std::optional<ReadStateSql> maybe_state_row = state_reader.next()){
                read_interupt_states = true;

                ReadStateSql state_row = maybe_state_row.value();
                temp_seed_state_map[state_row.seed][state_row.site_id] = state_row.degree_of_freedom;
            }


            SqlStatement<ReadCutoffSql> cutoff_statement(initial_state_database);
            SqlReader<ReadCutoffSql> cutoff_reader(cutoff_statement);
            while (std::optional<ReadCutoffSql> maybe_cutoff_row = cutoff_reader.next()){
                ReadCutoffSql cutoff_row = maybe_cutoff_row.value();

                temp_seed_step_map[cutoff_row.seed] = cutoff_row.step;
                temp_seed_time_map[cutoff_row.seed] = cutoff_row.time;
            }

            // bool read_trajectory_states = false;
            if (read_interupt_states == false) {
                // If the interupt_state table doesn't have entries, try to read from the trajectories table
                SqlStatement<ReadTrajectoriesSql> trajectory_statement(initial_state_database);
                SqlReader<ReadTrajectoriesSql> trajectory_reader(trajectory_statement);

                while (std::optional<ReadTrajectoriesSql> maybe_trajectory_row = trajectory_reader.next()) {
                    // read_trajectory_states = true;

                    ReadTrajectoriesSql trajectory_row = maybe_trajectory_row.value();
                    
                    Interaction* interaction = &model.all_interactions[trajectory_row.interaction_id];
                    temp_seed_state_map[trajectory_row.seed][trajectory_row.site_id_1] = interaction->right_state[0];
                    if (interaction->number_of_sites == 2) {
                        temp_seed_state_map[trajectory_row.seed][trajectory_row.site_id_2] = interaction->right_state[1];
                    }

                    if (trajectory_row.step > temp_seed_step_map[trajectory_row.seed]) {
                        temp_seed_step_map[trajectory_row.seed] = trajectory_row.step;
                        temp_seed_time_map[trajectory_row.seed] = trajectory_row.time;
                    }
                }
            }

            seed_state_map = temp_seed_state_map;
            seed_step_map = temp_seed_step_map;
            seed_time_map = temp_seed_time_map;

        };

    void run_dispatcher();
    void record_simulation_history(HistoryPacket history_packet);
    void record_state(StateHistoryPacket state_history_packet);
    void record_cutoff(CutoffHistoryPacket cutoff_history_packet);
};

void signalHandler(int signum) {
    
    write_error_message("received signal " + std::to_string(signum) + "\n");
    switch (signum) {
        case 2:
            // A SIGINT has been issued, handle this gracefully 
            // by writing current states and trajectories to the db

            write_error_message("SIGINT received. Terminating NPMC run(s) early.\n");
            do_shutdown = 1;
            shutdown_requested = true;

            write_error_message("Writing current states and history queue to the initial_states db\n");
            break;
        case 15:
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
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql>

void Dispatcher<Solver, Model, WriteTrajectoriesSql, ReadTrajectoriesSql, WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql>::run_dispatcher() {

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
                cutoff_history_queue,
                seed_queue,
                cutoff,
                running.begin() + i,
                seed_state_map,
                seed_step_map,
                seed_time_map
                )
            );
    }
    // Unset the sigmask so that the parent thread resumes catching errors as normal
    pthread_sigmask(SIG_UNBLOCK, &mask, NULL);
    
    action.sa_handler = signalHandler;
    sigemptyset(&action.sa_mask);
    action.sa_flags = 0;
    sigaction(SIGINT, &action, NULL);
    sigaction(SIGTERM, &action, NULL);

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
            // record_state(std::move(state_history_packet));
        };

        std::optional<CutoffHistoryPacket>
            maybe_cutoff_history_packet = cutoff_history_queue.get_history();

        if (maybe_cutoff_history_packet) {
            CutoffHistoryPacket cutoff_history_packet = std::move(maybe_cutoff_history_packet.value());
            record_cutoff(std::move(cutoff_history_packet));
        }
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
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql
    >
void Dispatcher<Solver, Model, WriteTrajectoriesSql, ReadTrajectoriesSql, WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql>::record_simulation_history(HistoryPacket history_packet) {
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
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql
    >
void Dispatcher<Solver, Model, WriteTrajectoriesSql, ReadTrajectoriesSql, WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql>::record_state(StateHistoryPacket state_history_packet) {

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

template <
    typename Solver,
    typename Model,
    typename WriteTrajectoriesSql,
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql
    >
void Dispatcher<Solver, Model, WriteTrajectoriesSql, ReadTrajectoriesSql, WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql>::record_cutoff(CutoffHistoryPacket cutoff_history_packet) {

    // Wipe the database of the states corresponding to this seed. This gets rid of the previously written states
    std::string delete_statement = "DELETE FROM interupt_cutoff WHERE seed = " + std::to_string(cutoff_history_packet.seed) + ";";
    
    initial_state_database.exec(delete_statement);
    
    initial_state_database.exec("BEGIN;");

    for (unsigned long int i = 0; i < cutoff_history_packet.cutoff.size(); i++) {
        cutoff_writer.insert(
            model.cutoff_history_element_to_sql(
                (int) cutoff_history_packet.seed,
                cutoff_history_packet.cutoff[i]));
    }
    initial_state_database.exec("COMMIT;");

    std::cerr << time_stamp()
              << "wrote cutoff for trajectory "
              << cutoff_history_packet.seed
              << " to database\n";
};
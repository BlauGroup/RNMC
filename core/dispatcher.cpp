#include <map>

#include "dispatcher.h"
#include "simulator_payload.h"

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename WriteTrajectoriesSql,
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql, 
    typename StateHistory, 
    typename TrajHistory, 
    typename CutoffHistory, 
    typename Sim, 
    typename State>

Dispatcher<Solver, Model, Parameters, WriteTrajectoriesSql,  ReadTrajectoriesSql, 
    WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql, StateHistory, TrajHistory,
    CutoffHistory, Sim, State>::Dispatcher(
        std::string model_database_file,
        std::string initial_state_database_file,
        unsigned long int number_of_simulations,
        unsigned long int base_seed,
        int number_of_threads,
        Cutoff cutoff,
        Parameters parameters):
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

            SeedQueue temp_seed_queue = SeedQueue(number_of_simulations, base_seed);
            std::map<int, State> temp_seed_state_map;
            std::map<int, int> temp_seed_step_map;
            std::map<int, double> temp_seed_time_map;

            SqlStatement<ReadStateSql> state_statement(initial_state_database);
            SqlReader<ReadStateSql> state_reader(state_statement);

            SqlStatement<ReadCutoffSql> cutoff_statement(initial_state_database);
            SqlReader<ReadCutoffSql> cutoff_reader(cutoff_statement);

            // If the interrupt_state table doesn't have entries, try to read from the trajectories table
            SqlStatement<ReadTrajectoriesSql> trajectory_statement(initial_state_database);
            SqlReader<ReadTrajectoriesSql> trajectory_reader(trajectory_statement);

            model.checkpoint(state_reader, cutoff_reader, trajectory_reader, 
                                                        temp_seed_state_map, temp_seed_step_map, 
                                                        temp_seed_queue, temp_seed_time_map, model);

           
            seed_state_map = temp_seed_state_map;
            seed_step_map = temp_seed_step_map;
            seed_time_map = temp_seed_time_map;

};


template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename WriteTrajectoriesSql,
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql, 
    typename StateHistory, 
    typename TrajHistory, 
    typename CutoffHistory, 
    typename Sim, 
    typename State>

void Dispatcher<Solver, Model, Parameters, WriteTrajectoriesSql,  ReadTrajectoriesSql, 
    WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql, StateHistory, TrajHistory,
    CutoffHistory, Sim, State>::signalHandler(int signum) {
    
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
    
} //signalHandler()

/* ---------------------------------------------------------------------- */

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename WriteTrajectoriesSql,
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql, 
    typename StateHistory, 
    typename TrajHistory, 
    typename CutoffHistory, 
    typename Sim, 
    typename State>

void Dispatcher<Solver, Model, Parameters, WriteTrajectoriesSql,  ReadTrajectoriesSql, 
    WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql, StateHistory, TrajHistory,
    CutoffHistory, Sim, State>::run_dispatcher() {
    
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
            [](SimulatorPayload<Solver, Model, StateHistory, TrajHistory, 
            CutoffHistory, Sim, State> payload) {payload.run_simulator();},
            SimulatorPayload<Solver, Model, StateHistory, TrajHistory, 
            CutoffHistory, Sim, State>(
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
    
    //action.sa_handler = signalHandler;
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


        std::optional<HistoryPacket<TrajHistory>>
            maybe_history_packet = history_queue.get_history();

        if (maybe_history_packet) {
            HistoryPacket<TrajHistory> history_packet = std::move(maybe_history_packet.value());
            record_simulation_history(std::move(history_packet));
        }

        std::optional<HistoryPacket<StateHistory>>
            maybe_state_history_packet = state_history_queue.get_history();

        if (maybe_state_history_packet) {
            HistoryPacket<StateHistory> state_history_packet = std::move(maybe_state_history_packet.value());
            record_state(std::move(state_history_packet));
        }

        std::optional<HistoryPacket<CutoffHistory>>
            maybe_cutoff_history_packet = cutoff_history_queue.get_history();

        if (maybe_cutoff_history_packet) {
            HistoryPacket<CutoffHistory> cutoff_history_packet = std::move(maybe_cutoff_history_packet.value());
            record_cutoff(std::move(cutoff_history_packet));
        }
    }

    for (int i = 0; i < number_of_threads; i++) threads[i].join();

    initial_state_database.exec(
        "DELETE FROM trajectories WHERE rowid NOT IN"
        "(SELECT MIN(rowid) FROM trajectories GROUP BY seed, step);");

    std::cerr << time::time_stamp()
              << "removing duplicate trajectories...\n";

} //run_dispatcher()

/* ---------------------------------------------------------------------- */

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename WriteTrajectoriesSql,
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql, 
    typename StateHistory, 
    typename TrajHistory, 
    typename CutoffHistory, 
    typename Sim, 
    typename State>

void Dispatcher<Solver, Model, Parameters, WriteTrajectoriesSql,  ReadTrajectoriesSql, 
    WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql, StateHistory, TrajHistory,
    CutoffHistory, Sim, State>::record_simulation_history(HistoryPacket<TrajHistory> history_packet) {
    initial_state_database.exec("BEGIN;");

    for (unsigned long int i = 0; i < history_packet.history.size(); i++) {
        trajectories_writer.insert(
            model.history_element_to_sql(
                (int) history_packet.seed,
                history_packet.history[i]));


    }
    initial_state_database.exec("COMMIT;");

    std::cerr << time::time_stamp()
              << "wrote "
              << history_packet.history.size()
              << " events from trajectory "
              << history_packet.seed
              << " to database\n";

} //record_simulation_history()

/* ---------------------------------------------------------------------- */

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename WriteTrajectoriesSql,
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql, 
    typename StateHistory, 
    typename TrajHistory, 
    typename CutoffHistory, 
    typename Sim, 
    typename State>

void Dispatcher<Solver, Model, Parameters, WriteTrajectoriesSql,  ReadTrajectoriesSql, 
    WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql, StateHistory, TrajHistory,
    CutoffHistory, Sim, State>::record_state(HistoryPacket<StateHistory> state_history_packet) {

    // Wipe the database of the states corresponding to this seed. This gets rid of the previously written states
    std::string delete_statement = "DELETE FROM interrupt_state WHERE seed = " + std::to_string(state_history_packet.seed) + ";";
    
    initial_state_database.exec(delete_statement);
    
    initial_state_database.exec("BEGIN;");

    for (unsigned long int i = 0; i < state_history_packet.history.size(); i++) {
        state_writer.insert(
            model.state_history_element_to_sql(
                (int) state_history_packet.seed,
                state_history_packet.history[i]));
    }
    initial_state_database.exec("COMMIT;");

    std::cerr << time::time_stamp()
              << "wrote "
              << state_history_packet.history.size()
              << " states for trajectory "
              << state_history_packet.seed
              << " to database\n";
} //record_state()

/* ---------------------------------------------------------------------- */

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename WriteTrajectoriesSql,
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql, 
    typename StateHistory, 
    typename TrajHistory, 
    typename CutoffHistory, 
    typename Sim, 
    typename State>

void Dispatcher<Solver, Model, Parameters, WriteTrajectoriesSql,  ReadTrajectoriesSql, 
    WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql, StateHistory, TrajHistory,
    CutoffHistory, Sim, State>::record_cutoff(HistoryPacket<CutoffHistory> cutoff_history_packet) {

    // Wipe the database of the states corresponding to this seed. This gets rid of the previously written states
    std::string delete_statement = "DELETE FROM interrupt_cutoff WHERE seed = " + std::to_string(cutoff_history_packet.seed) + ";";
    
    initial_state_database.exec(delete_statement);
    
    initial_state_database.exec("BEGIN;");

    for (unsigned long int i = 0; i < cutoff_history_packet.history.size(); i++) {
        cutoff_writer.insert(
            model.cutoff_history_element_to_sql(
                (int) cutoff_history_packet.seed,
                cutoff_history_packet.history[i]));
    }
    initial_state_database.exec("COMMIT;");

    std::cerr << time::time_stamp()
              << "wrote cutoff for trajectory "
              << cutoff_history_packet.seed
              << " to database\n";
} //record_cutoff()

/* ---------------------------------------------------------------------- */

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename WriteTrajectoriesSql,
    typename ReadTrajectoriesSql,
    typename WriteStateSql,
    typename ReadStateSql,
    typename WriteCutoffSql,
    typename ReadCutoffSql, 
    typename StateHistory, 
    typename TrajHistory, 
    typename CutoffHistory, 
    typename Sim, 
    typename State>

void Dispatcher<Solver, Model, Parameters, WriteTrajectoriesSql,  ReadTrajectoriesSql, 
    WriteStateSql, ReadStateSql, WriteCutoffSql, ReadCutoffSql, StateHistory, TrajHistory,
    CutoffHistory, Sim, State>::write_error_message(std::string s){
    char char_array[s.length()+1];
    strcpy(char_array, s.c_str());

    write(STDERR_FILENO, char_array, sizeof(char_array) - 1);

} // write_error_message()
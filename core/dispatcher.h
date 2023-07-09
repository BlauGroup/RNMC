#ifndef RNMC_DISPATCHER_H
#define RNMC_DISPATCHER_H

#include <mutex>
#include <thread>
#include <csignal>
#include <iostream>
#include <string>
#include <map>

#include "sql.h"
#include "queues.h"
#include "RNMC_types.h"
#include "simulation.h"
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

class Dispatcher {

public:
    SqlConnection model_database;
    SqlConnection initial_state_database;
    Model model;
    SqlStatement<WriteTrajectoriesSql> trajectories_stmt;
    SqlWriter<WriteTrajectoriesSql> trajectories_writer;
    
    SqlStatement<WriteStateSql> state_stmt;
    SqlWriter<WriteStateSql> state_writer;
    
    SqlStatement<WriteCutoffSql> cutoff_stmt;
    SqlWriter<WriteCutoffSql> cutoff_writer;

    HistoryQueue<HistoryPacket<TrajHistory>> history_queue;
    HistoryQueue<HistoryPacket<StateHistory>> state_history_queue;
    HistoryQueue<HistoryPacket<CutoffHistory>> cutoff_history_queue;

    SeedQueue seed_queue;
    std::vector<std::thread> threads;
    std::vector<bool> running;
    Cutoff cutoff;
    TypeOfCutoff type_of_cutoff;
    int number_of_simulations;
    int number_of_threads;
    struct sigaction action;
    std::map<int, State> seed_state_map;
    std::map<int, int> seed_step_map;
    std::map<int, double> seed_time_map;
    
    
    Dispatcher(
        std::string model_database_file,
        std::string initial_state_database_file,
        unsigned long int number_of_simulations,
        unsigned long int base_seed,
        int number_of_threads,
        Cutoff cutoff,
        Parameters parameters);

    void static signalHandler(int signum);
    void run_dispatcher();
    void record_simulation_history(HistoryPacket<TrajHistory> traj_history_packet);
    void record_state(HistoryPacket<StateHistory> state_history_packet);
    void record_cutoff(HistoryPacket<CutoffHistory> cutoff_history_packet);
    void static write_error_message(std::string s);
    
};

#include "dispatcher.cpp"

#endif 
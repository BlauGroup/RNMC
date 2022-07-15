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

constexpr int history_chunk_size = 20000;

template <typename Solver, typename Model>
struct SimulatorPayload {
    Model &model;
    HistoryQueue<HistoryPacket> &history_queue;
    SeedQueue &seed_queue;
    Cutoff cutoff;
    std::vector<bool>::iterator running;
    bool isLGMC;

    SimulatorPayload(
        Model &model,
        HistoryQueue<HistoryPacket> &history_queue,
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


            Simulation<Solver, Model> simulation (model, seed, history_chunk_size, history_queue);


            switch(cutoff.type_of_cutoff) {
            case step_termination :
                simulation.execute_steps(cutoff.bound.step, isLGMC);
                break;
            case time_termination :
                simulation.execute_time(cutoff.bound.time, isLGMC);
                break;
            }


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
    std::vector<bool> running;
    Cutoff cutoff;
    TypeOfCutoff type_of_cutoff;
    int number_of_simulations;
    int number_of_threads;
    bool isLGMC;

    Dispatcher(std::string file_in); 
    void run_dispatcher();
    void record_simulation_history(HistoryPacket history_packet);
};


template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename TrajectoriesSql>

Dispatcher::Dispatcher<Solver, Model, Parameters, TrajectoriesSql>Dispatcher(std::string file_in, bool isLGMC) {

    std::string reaction_database;
    std::string initial_state_database;
    unsigned long int number_of_simulations;
    unsigned long int number_of_threads;
    int base_seed;
    int cutoff;
    char cutoff_flag;
    bool isLGMC_in;

    Cutoff cutoff = {
        .bound =  { .step =  0 },
        .type_of_cutoff = step_termination
    };

    std::cin >> reaction_database >> initial_state_database
    >> number_of_simulations >> base_seed >> number_of_threads >> 
    cutoff_flag >> cutoff >> isLGMC_in;

    if(std::cin.fail()) {
        std::cout << "Incorrect file arguments.\n";
        exit(EXIT_FAILURE);
    }

    if(cutoff_flag == 'S') {
        cutoff.bound.step = cutoff;
        cutoff.type_of_cutoff = step_termination;
        break;
    }
    else if(cutoff_flag == 'T') {
        cutoff.bound.time = atof(optarg);
        cutoff.type_of_cutoff = time_termination;
        break;
    }
    else {
        std::cout << "Incorrect cutoff flag (S/T).";
        exit(EXIT_FAILURE);
    }

    model_database = SQLConnection (model_database, SQLITE_OPEN_READWRITE);
    initial_state_database = SqlConnection (initial_state_database, SQLITE_OPEN_READWRITE);


    trajectories_stmt = SqlStatement<TrajectoriesSql> (initial_state_database);
    trajectories_writer = SqlWriter<TrajectoriesSql> (trajectories_stmt);
    history_queue = HistoryQueue<HistoryPacket> ();
    seed_queue = SeedQueue (number_of_simulations, base_seed);

    running = std::vector<bool> (number_of_threads, false);
    cutoff (cutoff),
    number_of_simulations = number_of_simulations;
    number_of_threads = number_of_threads;
    isLGMC = isLGMC_in;

    if(isLGMC) {
        float latconst;                               
        int boxxlo,boxxhi,boxylo,                   
        boxyhi,boxzlo,boxzhi;                       
        float xperiodic,yperiodic,zperiodic;

        std::cin >> latconst >> boxxlo >> boxxhi >> boxylo >> boxyhi >> boxzlo >> boxzhi
        >> xperiodic >> yperiodic >> zperiodic;

        if(std::cin.fail()) {
            std::cout << "Incorrect file arguments.\n";
            exit(EXIT_FAILURE);
        }   

        lattice_reaction_network = LatticeReactionNetwork(model_database, initial_state_database, 
        parameters);

        // create lattice
        lattice = new Lattice(latconst, boxxlo, boxxhi, boxylo,
                          boxyhi, boxzlo, boxzhi, xperiodic, yperiodic, zperiodic);

        model = new LGMC(lattice, lattice_reaction_network);

    } else {
        model = Model (model_database, initial_state_database, parameters);
    }

    model.init(model_database);

}

template <
    typename Solver,
    typename Model,
    typename Parameters,
    typename TrajectoriesSql>

void Dispatcher<Solver, Model, Parameters, TrajectoriesSql>::run_dispatcher() {


    threads.resize(number_of_threads);
    for (int i = 0; i < number_of_threads; i++) {
        running[i] = true;
        threads[i] = std::thread (
            [](SimulatorPayload<Solver, Model> payload) {payload.run_simulator();},
            SimulatorPayload<Solver, Model> (
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


        std::optional<HistoryPacket>
            maybe_history_packet = history_queue.get_history();

        if (maybe_history_packet) {
            HistoryPacket history_packet = std::move(maybe_history_packet.value());
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
    typename TrajectoriesSql
    >
void Dispatcher<Solver, Model, Parameters, TrajectoriesSql>::record_simulation_history(HistoryPacket history_packet) {
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

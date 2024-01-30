#include <getopt.h>

#include "../core/nano_particle_simulation.h"
#include "../core/dispatcher.h"
#include "nano_particle.h"

void print_usage() {
    std::cout << "Usage: specify the following options\n"
              << "--nano_particle_database\n"
              << "--initial_state_database\n"
              << "--number_of_simulations\n"
              << "--base_seed\n"
              << "--thread_count\n"
              << "--step_cutoff|time_cutoff\n"
              << "--checkpoint\n";

} // print_usage()

/* ---------------------------------------------------------------------- */

int main(int argc, char **argv) {
    if (argc != 8) {
        print_usage();
        exit(EXIT_FAILURE);
    }

    struct option long_options[] = {
        {"nano_particle_database", required_argument, NULL, 1},
        {"initial_state_database", required_argument, NULL, 2},
        {"number_of_simulations", required_argument, NULL, 3},
        {"base_seed", required_argument, NULL, 4},
        {"thread_count", required_argument, NULL, 5},
        {"step_cutoff", optional_argument, NULL, 6},
        {"time_cutoff", optional_argument, NULL, 7},
        {"checkpoint", required_argument, NULL, 8},
        {NULL, 0, NULL, 0}
        // last element of options array needs to be filled with zeros
    };

    int c;
    int option_index = 0;

    char *nano_particle_database = nullptr;
    char *initial_state_database = nullptr;
    int number_of_simulations = 0;
    int base_seed = 0;
    int thread_count = 0;
    bool isCheckpoint = false;
    
    Cutoff cutoff = {
        .bound =  { .step =  0 },
        .type_of_cutoff = step_termination
    };

    while ((c = getopt_long_only(
                argc, argv, "",
                long_options,
                &option_index)) != -1) {

        switch (c) {

        case 1:
            nano_particle_database = optarg;
            break;

        case 2:
            initial_state_database = optarg;
            break;

        case 3:
            number_of_simulations = atoi(optarg);
            break;

        case 4:
            base_seed = atoi(optarg);
            break;

        case 5:
            thread_count = atoi(optarg);
            break;

        case 6:
            cutoff.bound.step = atoi(optarg);
            cutoff.type_of_cutoff = step_termination;
            break;

        case 7:
            cutoff.bound.time = atof(optarg);
            cutoff.type_of_cutoff = time_termination;
            break;

        case 8:
            isCheckpoint = atof(optarg);
            break;

        default:
            // if an unexpected argument is passed, exit
            print_usage();
            exit(EXIT_FAILURE);
            break;

        }

    }

    NanoParticleParameters parameters {
        .isCheckpoint = isCheckpoint
    };
    
    Dispatcher<
    NanoSolver,
    NanoParticle,
    NanoParticleParameters,
    NanoWriteTrajectoriesSql,
    NanoReadTrajectoriesSql,
    NanoWriteStateSql,
    NanoReadStateSql,
    WriteCutoffSql,
    ReadCutoffSql, 
    NanoStateHistoryElement, 
    NanoTrajectoryHistoryElement, 
    CutoffHistoryElement,
    NanoParticleSimulation, 
    std::vector<int>>

        dispatcher (
            nano_particle_database,
            initial_state_database,
            number_of_simulations,
            base_seed,
            thread_count,
            cutoff, 
            parameters
            );

    dispatcher.run_dispatcher();
    exit(EXIT_SUCCESS);

}
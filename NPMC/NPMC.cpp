#include <getopt.h>
#include "../core/dispatcher.h"
#include "sql_types.h"
#include "nano_particle.h"

void print_usage() {
    std::cout << "Usage: specify the following options\n"
              << "--nano_particle_database\n"
              << "--initial_state_database\n"
              << "--number_of_simulations\n"
              << "--base_seed\n"
              << "--thread_count\n"
              << "--step_cutoff\n";
}

int main(int argc, char **argv) {
    if (argc != 7) {
        print_usage();
        exit(EXIT_FAILURE);
    }


    struct option long_options[] = {
        {"nano_particle_database", required_argument, NULL, 1},
        {"initial_state_database", required_argument, NULL, 2},
        {"number_of_simulations", required_argument, NULL, 3},
        {"base_seed", required_argument, NULL, 4},
        {"thread_count", required_argument, NULL, 5},
        {"step_cutoff", required_argument, NULL, 6},
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
    int step_cutoff = 0;

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
            step_cutoff = atoi(optarg);
            break;

        default:
            // if an unexpected argument is passed, exit
            print_usage();
            exit(EXIT_FAILURE);
            break;

        }

    }
    NanoParticleParameters parameters = {};

    Dispatcher<
        LinearSolver,
        NanoParticle,
        NanoParticleParameters,
        TrajectoriesSql
        >

        dispatcher (
            nano_particle_database,
            initial_state_database,
            number_of_simulations,
            base_seed,
            thread_count,
            step_cutoff,
            parameters
            );

    dispatcher.run_dispatcher();
    exit(EXIT_SUCCESS);

}

#include <getopt.h>

#include "sql_types.h"
#include "linear_solver.h"
#include "../core/dispatcher.h"
#include "../core/reaction_network_simulation.h"
#include "../core/energy_reaction_network_simulation.h"

void print_usage() {
    std::cout << "Usage: specify the following options\n"
              << "--reaction_database\n"
              << "--initial_state_database\n"
              << "--number_of_simulations\n"
              << "--base_seed\n"
              << "--thread_count\n"
              << "--step_cutoff|time_cutoff\n"
              << "--energy_budget\n"
              << "--checkpoint\n";
} // print_usage()

/*---------------------------------------------------------------------------*/

int main(int argc, char **argv) {
    // Here we just check if there are 7 or 8 args.
    // This is a lazy way to let energy budget be optional.
    // This could potentially break if both step and time cutoff are specified while energy budget is not
    if (argc < 8 || argc > 9) {
        print_usage();
        exit(EXIT_FAILURE);
    }

    struct option long_options[] = {
        {"reaction_database", required_argument, NULL, 1},
        {"initial_state_database", required_argument, NULL, 2},
        {"number_of_simulations", required_argument, NULL, 3},
        {"base_seed", required_argument, NULL, 4},
        {"thread_count", required_argument, NULL, 5},
        {"step_cutoff", optional_argument, NULL, 6},
        {"time_cutoff", optional_argument, NULL, 7},
        {"energy_budget", optional_argument, NULL, 8},
        {"checkpoint", required_argument, NULL, 9},
        {NULL, 0, NULL, 0}
        // last element of options array needs to be filled with zeros
    };

    int c;
    int option_index = 0;

    char *reaction_database = nullptr;
    char *initial_state_database = nullptr;
    int number_of_simulations = 0;
    int base_seed = 0;
    int thread_count = 0;
    double energy_budget = 0;
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
            reaction_database = optarg;
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
            energy_budget = atof(optarg);
            break;

        case 9:
            isCheckpoint = atof(optarg);
            std::cout << isCheckpoint << std::endl;
            break;

        default:
            // if an unexpected argument is passed, exit
            print_usage();
            exit(EXIT_FAILURE);
            break;
        }
    }

    if(energy_budget == 0) {
        ReactionNetworkParameters parameters {
            .isCheckpoint = isCheckpoint
        };

        Dispatcher<
            LinearSolver,
            GillespieReactionNetwork,
            ReactionNetworkParameters,
            ReactionNetworkWriteTrajectoriesSql,
            ReactionNetworkReadTrajectoriesSql,
            ReactionNetworkWriteStateSql,
            ReactionNetworkReadStateSql,
            WriteCutoffSql,
            ReadCutoffSql, 
            ReactionNetworkStateHistoryElement, 
            ReactionNetworkTrajectoryHistoryElement, 
            CutoffHistoryElement, 
            ReactionNetworkSimulation<LinearSolver>, 
            std::vector<int>>

        dispatcher (
        reaction_database,
        initial_state_database,
        number_of_simulations,
        base_seed,
        thread_count,
        cutoff,
        parameters
        );

        dispatcher.run_dispatcher();
        
    }
    else {
        EnergyReactionNetworkParameters parameters{
            .energy_budget = energy_budget,
            .isCheckpoint = isCheckpoint
        };

        Dispatcher<
            LinearSolver,
            EnergyReactionNetwork,
            EnergyReactionNetworkParameters,
            ReactionNetworkWriteTrajectoriesSql,
            ReactionNetworkReadTrajectoriesSql,
            ReactionNetworkWriteStateSql,
            ReactionNetworkReadStateSql,
            EnergyNetworkWriteCutoffSql,
            EnergyNetworkReadCutoffSql, 
            ReactionNetworkStateHistoryElement, 
            ReactionNetworkTrajectoryHistoryElement, 
            EnergyNetworkCutoffHistoryElement, 
            EnergyReactionNetworkSimulation<TreeSolver>, 
            EnergyState>

        dispatcher (
        reaction_database,
        initial_state_database,
        number_of_simulations,
        base_seed,
        thread_count,
        cutoff,
        parameters
        );

        dispatcher.run_dispatcher();

    }

    exit(EXIT_SUCCESS);
}
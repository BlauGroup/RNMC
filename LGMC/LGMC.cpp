#include "../core/dispatcher.h"
#include "../GMC/sql_types.h"
#include "lattice_reaction_network.h"

#include <string>
#include <fstream>
#include <iostream>
#include <getopt.h>

void print_usage() {
    
    std::cout << "Usage: specify the following in the input file\n"
              << "reaction_database\n"
              << "initial_state_database\n"
              << "number_of_simulations\n"
              << "base_seed\n"
              << "step_cutoff\n"
              << "parameters\n";
    
} // print_usage()

void print_usage_LGMC_parameters() {
 
     std::cout << "Usage: specify the following in the input file\n"
              << "lattice constant\n"
              << "box x lower boundary\n"
              << "box x upper boundary\n"
              << "box y lower boundary\n"
              << "box y upper boundary\n"
              << "box z lower boundary\n"
              << "box z upper boundary\n"
              << "Periodicity in x dimension (T/F)\n"
              << "Periodicity in y dimension (T/F)\n"
              << "Periodicity in z dimension (T/F)\n"
              << "Temperature\n"
              << "Electron free energy\n"
              << "Charge transfer style (M/B)\n"
              << "Is Add Site (T/F)\n";
}

/* ---------------------------------------------------------------------- */

int main(int argc, char **argv) {
    if (argc != 7) {
        print_usage();
        exit(EXIT_FAILURE);
    }


    struct option long_options[] = {
        {"reaction_database", required_argument, NULL, 1},
        {"initial_state_database", required_argument, NULL, 2},
        {"number_of_simulations", required_argument, NULL, 3},
        {"base_seed", required_argument, NULL, 4},
        //{"thread_count", required_argument, NULL, 5},
        {"step_cutoff", optional_argument, NULL, 6},
        {"time_cutoff", optional_argument, NULL, 7},
        {"parameters", required_argument, NULL, 8},
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
    std::string LGMC_params_file; 

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
            LGMC_params_file = optarg;
            break;

        default:
            // if an unexpected argument is passed, exit
            print_usage();
            exit(EXIT_FAILURE);
            break;

        }

    }                    

    // read in LGMC parameters from file 
    std::ifstream fin;
    fin.open(LGMC_params_file);

    if(!fin.is_open()) {
        std::cout << "Failed to open file: " << LGMC_params_file << "\n";
        print_usage_LGMC_parameters();
        exit(EXIT_FAILURE);
    }

    float latconst;                               
    float boxxlo,boxxhi,boxylo,                   
          boxyhi,boxzlo,boxzhi;                       
    char xperiod, yperiod, zperiod;
    bool xperiodic, yperiodic, zperiodic;
    float temperature;
    float g_e;
    bool is_add_site;
    char add_site;
    ChargeTransferStyle charge_transfer_style;
    char ct_style;


    std::cin >> latconst >> boxxlo >> boxxhi >> boxylo >> boxyhi >> boxzlo >> boxzhi
    >> xperiod >> yperiod >> zperiod >> temperature >> g_e >> add_site >> ct_style;

    if(std::cin.fail()) {
        std::cout << "Incorrect parameter file arguments.\n";
        exit(EXIT_FAILURE);
    }

    if(add_site == 'T') {
        is_add_site = true;
    }
    else if (xperiod == 'F') {
        is_add_site = false;
    }else {
        std::cout << "Incorrect parameter file argument for add site.\n";
        exit(EXIT_FAILURE);
    }

    if(xperiod == 'T') {
        xperiodic = true;
    }
    else if (xperiod == 'F') {
        xperiodic = false;
    }else {
        std::cout << "Incorrect parameter file argument for x periodicity.\n";
        exit(EXIT_FAILURE);
    }

    if(yperiod == 'T') {
        yperiodic = true;
    }
    else if (yperiod == 'F') {
        yperiodic = false;
    }else {
        std::cout << "Incorrect parameter file argument for y periodicity.\n";
        exit(EXIT_FAILURE);
    }

    if(zperiod == 'T') {
        zperiodic = true;
    }
    else if (zperiod == 'F') {
        zperiodic = false;
    }else {
        std::cout << "Incorrect parameter file argument for z periodicity.\n";
        exit(EXIT_FAILURE);
    }

    if(ct_style == 'M') {
        charge_transfer_style = ChargeTransferStyle::MARCUS;
    }
    else if (ct_style == 'B') {
        charge_transfer_style = ChargeTransferStyle::BUTLER_VOLMER;
    }else {
        std::cout << "Incorrect parameter file argument for charge transfer style.\n";
        exit(EXIT_FAILURE);
    }


    LGMCParameters parameters{.latconst = latconst, .boxxlo = boxxlo, .boxxhi = boxxhi, 
                              .boxylo = boxyhi, .boxzlo = boxzlo, .boxzhi = boxzhi, 
                              .xperiodic = xperiodic, .yperiodic = yperiodic, .zperiodic = zperiodic, 
                              .temperature = temperature, .g_e = g_e, .is_add_sites = is_add_site,
                              .charge_transfer_style = charge_transfer_style};                               

    Dispatcher<
        LatSolver,
        LGMC_NS::LatticeReactionNetwork,
        LGMCParameters,
        TrajectoriesSql
        >

        dispatcher (
        reaction_database,
        initial_state_database,
        number_of_simulations,
        base_seed,
        thread_count,
        cutoff,
        parameters, 
        true
        );

    dispatcher.run_dispatcher();
    exit(EXIT_SUCCESS);


}
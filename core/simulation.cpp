#include <string>
#include "simulation.h"

template <typename Solver>
void Simulation<Solver>::execute_steps(int step_cutoff) {
    while(execute_step()) {
        if (this->step > step_cutoff) {
            break;
        } else if (do_shutdown || shutdown_requested.load()) {
            // Handle shutdown request from SIGTERM
            write_error_message("Received termination request on thread- cleaning up\n");
            break;
        }
    }
} // execute_steps()

/* ---------------------------------------------------------------------- */

template <typename Solver>
void Simulation<Solver>::execute_time(double time_cutoff) {
    while(execute_step()) {
        if (time > time_cutoff) {
            break;
        } else if (do_shutdown || shutdown_requested.load()) {
            // Handle shutdown request from SIGTERM
            write_error_message("Received termination request on thread- cleaning up\n");
            break;
        }
    }
} // execute_time()

/* ---------------------------------------------------------------------- */

template <typename Solver>
bool Simulation<Solver>::execute_step() {
    return false;
} // execute_steps()

/* ---------------------------------------------------------------------- */

template <typename Solver>
void Simulation<Solver>::write_error_message(std::string s){
    char char_array[s.length()+1];
    strcpy(char_array, s.c_str());

    write(STDERR_FILENO, char_array, sizeof(char_array) - 1);
} // write_error_message()
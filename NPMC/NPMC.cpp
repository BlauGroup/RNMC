#include "../core/sql.h"
#include "../core/solvers.h"
#include "simulation.h"
#include <functional>

int main() {


    SqlConnection nanoparticle_database (
        "./test_materials/NPMC/np.sqlite",
        SQLITE_OPEN_READONLY);

    SqlConnection initial_state_database (
        "./test_materials/NPMC/initial_state.sqlite",
        SQLITE_OPEN_READWRITE);

    NanoParticle nano_particle (
        std::ref(nanoparticle_database),
        std::ref(initial_state_database));


    Simulation<TreeSolver> simulation (
        std::ref(nano_particle),
        42,
        1000);

    simulation.execute_steps(1000);
}

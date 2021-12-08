#include "../core/sql.h"
#include "nano_particle.h"
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
}

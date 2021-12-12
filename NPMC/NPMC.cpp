#include "../core/dispatcher.h"
#include "sql_types.h"
#include "nano_particle.h"


int main() {

    NanoParticleParameters parameters = {};

    Dispatcher<
        TreeSolver,
        NanoParticle<TreeSolver>,
        NanoParticleParameters,
        TrajectoriesSql
        >

        dispatcher (
            "./test_materials/NPMC/np.sqlite",
            "./test_materials/NPMC/initial_state.sqlite",
            1000,
            1000,
            8,
            2000,
            parameters
            );

    dispatcher.run_dispatcher();
    exit(EXIT_SUCCESS);

}

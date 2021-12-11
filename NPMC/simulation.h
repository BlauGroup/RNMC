#include "nano_particle.h"


struct HistoryElement : Reaction {
    double time;
};

template <typename Solver>
struct Simulation {
    NanoParticle &nano_particle;
    unsigned long int seed;
    std::vector<int> state;
    double time;
    int step;       // number of reactions which have occoured
    Solver solver;
    std::vector<HistoryElement> history;


    Simulation( NanoParticle &nano_particle,
                unsigned long int seed,
                int step_cutoff
        ) :
        nano_particle (nano_particle),
        seed (seed),
        state (nano_particle.initial_state),
        time (0.0),
        step (0),
        solver (seed, std::ref(nano_particle.initial_propensities)),
        history (step_cutoff + 1)
        {};
};


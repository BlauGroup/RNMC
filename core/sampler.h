#ifndef RNMC_SAMPLER_H
#define RNMC_SAMPLER_H

#include <gsl/gsl_rng.h>
#include <utility>

// we are using GSL random number generation because i don't trust
// random number generation to be consistent across various C++ stdlib
// implementations, and for the kind of MC simulator we are writing here,
// we want to be able to run it deterministically for testing purposes.

class Sampler {
private:
    gsl_rng *internal_rng_state;
public:
    unsigned long int seed;
    double generate() {
        return gsl_rng_uniform_pos(internal_rng_state);
    };

    Sampler(unsigned long int n) :
        seed (n) {
        internal_rng_state = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(internal_rng_state, seed);
        };

    ~Sampler() {
        gsl_rng_free(internal_rng_state);
    };

    // since we don't have access to gsl internal state, can't write
    // copy constructor
    Sampler(Sampler &other) = delete;

    // move constructor
    Sampler(Sampler &&other) :
        internal_rng_state (std::exchange(other.internal_rng_state, nullptr)),
        seed (other.seed)
        {};

    // since we don't have access to gsl internal state, can't write
    // copy assignment operator
    Sampler &operator=(Sampler &other) = delete;

    // move assignment operator
    Sampler &operator=(Sampler &&other) {
        seed = other.seed;

        // we move the existing internal_rng_state into other so
        // it gets freed when other is dropped.
        std::swap(internal_rng_state, other.internal_rng_state);
        return *this;
    };
};

#endif
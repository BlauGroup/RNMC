#include <gsl/gsl_rng.h>
#include <utility>
#include <iostream>

class Sampler {
private:
    gsl_rng *internal_rng_state;
public:
    unsigned long int seed;
    double generate() {
        return gsl_rng_uniform_pos(internal_rng_state);
    };

    Sampler(unsigned long int n) :
        seed{n} {
        internal_rng_state = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(internal_rng_state, seed);
        };

    ~Sampler() {
        gsl_rng_free(internal_rng_state);
    };

    // since we don't have access to gsl internal state, can't write
    // copy constructor
    Sampler(Sampler &other) = delete;

    Sampler(Sampler &&other) :
        seed{other.seed},
        internal_rng_state{std::exchange(other.internal_rng_state, nullptr)} {};

    // since we don't have access to gsl internal state, can't write
    // copy assignment operator
    Sampler &operator=(Sampler &other) = delete;

    Sampler &operator=(Sampler &&other) {
        seed = other.seed;
        std::swap(internal_rng_state, other.internal_rng_state);
        return *this;
    };

};


int main() {

    Sampler sampler(42);

    std::cout << sampler.generate() << sampler.generate();
}

#include "sampler.h"
#include <iostream>

int main() {

    Sampler sampler(42);
    int i;
    for (i = 0; i < 1000; i++) {
        std::cout << sampler.generate() << '\n';
    }
}


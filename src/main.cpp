#include "sampler.h"
#include <iostream>

int main() {

    Sampler sampler(42);
    Sampler foo(std::move(sampler));
    int i;
    for (i = 0; i < 1000; i++) {
        std::cout << foo.generate() << '\n';
    }
}


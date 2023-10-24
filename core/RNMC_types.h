#ifndef RNMC_RNMC_TYPES_H
#define RNMC_RNMC_TYPES_H

#include <vector>

enum TypeOfCutoff {step_termination, time_termination };

struct Cutoff {
    union  { int step; double time; } bound;
    TypeOfCutoff type_of_cutoff;
};

template <typename T>
struct HistoryPacket {
    unsigned long int seed;
    std::vector<T> history;
};

struct Update {
    unsigned long int index;
    double propensity;
};

struct Event {
    unsigned long int index;
    double dt;
};

#endif
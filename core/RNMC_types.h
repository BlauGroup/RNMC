#ifndef RNMC_TYPES_H
#define RNMC_TYPES_H

enum TypeOfCutoff {step_termination, time_termination };

struct Cutoff {
    union  { int step; double time; } bound;
    TypeOfCutoff type_of_cutoff;
};

template <typename T>
struct HistoryPacket {
    std::vector<T> history;
    unsigned long int seed;
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
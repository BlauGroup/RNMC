#pragma once
#include <vector>

struct Update {
    int index;
    double propensity;
};

struct Event {
    int index;
    double dt;
};

class Solver {
    virtual void update(Update update);
    virtual void update(std::vector<Update> updates);
    virtual Event event();
    virtual double get_propensity(int index);
    virtual double get_propensity_sum();
};


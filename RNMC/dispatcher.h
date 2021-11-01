#pragma once
#include <queue>
#include <mutex>
#include "simulation.h"


struct SeedQueue {
    std::queue<unsigned long int> seeds;
    std::mutex mutex;

    SeedQueue(unsigned long int number_of_seeds, unsigned long int base_seed) {
        for (unsigned long int i = base_seed; i < number_of_seeds; i++) {
            seeds.push(i);
        }
    }

    std::optional<unsigned long int> get_seed() {
        std::lock_guard<std::mutex> lock (mutex);

        if (seeds.empty()) {
            return std::optional<unsigned long int> ();
        } else {
            unsigned long int result = seeds.front();
            seeds.pop();
            return std::optional<unsigned long int> (result);
        }
    }
};

struct HistoryQueue {
    std::queue<std::vector<HistoryElement>> histories;
    std::mutex mutex;

    void insert_history(std::vector<HistoryElement> &&history) {
        std::lock_guard<std::mutex> lock (mutex);
        histories.push(std::move(history));
    }

    std::vector<HistoryElement> get_history() {
        std::lock_guard<std::mutex> lock (mutex);
        std::vector<HistoryElement> result = std::move(histories.front());
        histories.pop();
        return result;
    };

};


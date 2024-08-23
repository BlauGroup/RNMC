/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#ifndef RNMC_QUEUES_H
#define RNMC_QUEUES_H

#include <queue>
#include <mutex>
#include <optional>

struct SeedQueue
{
    std::queue<unsigned long int> seeds;
    std::mutex mutex;

    SeedQueue(unsigned long int number_of_seeds, unsigned long int base_seed)
    {
        for (unsigned long int i = base_seed;
             i < number_of_seeds + base_seed;
             i++)
        {
            seeds.push(i);
        }
    }

    std::optional<unsigned long int> get_seed()
    {
        std::lock_guard<std::mutex> lock(mutex);

        if (seeds.empty())
        {
            return std::optional<unsigned long int>();
        }
        else
        {
            unsigned long int result = seeds.front();
            seeds.pop();
            return std::optional<unsigned long int>(result);
        }
    }
};

/* ------------------------------------------------------------------- */

template <typename T>
struct HistoryQueue
{
    // the flow of trajectory histories from the simulator threads to
    // the dispatcher is subtle and important. The vector of histories
    // is allocated by the simulator. Once the simulation is finished,
    // it is moved into the history queue. Then the dispatcher moves
    // it out of the history queue, writes it into the initial state
    // database and then frees it. Without all the carefully placed
    // std::move calls, there will be a lot of unnecessary allocations
    // and frees. In typical C++ fashion, this is made significantly
    // more complex by the subtle move semantics of std::optional. If
    // you are going to change this code, you need to use gdb to check
    // that you haven't accidently introduced extra allocations and
    // frees (i.e the vector allocated by the simulation thread points
    // to exactly the same memory as the vector which is used to write
    // to the initial state database, and every time a move is
    // supposed to happen, the old reference is actually zerod out).
    std::queue<T> history_packets;
    std::mutex mutex;

    bool empty()
    {
        std::lock_guard<std::mutex> lock(mutex);
        return history_packets.empty();
    }

    void insert_history(T history_packet)
    {
        std::lock_guard<std::mutex> lock(mutex);
        history_packets.push(std::move(history_packet));
    }

    std::optional<T> get_history()
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (history_packets.empty())
        {
            return std::optional<T>();
        }
        else
        {
            T result = std::move(history_packets.front());
            history_packets.pop();
            return std::optional<T>(std::move(result));
        }
    };
};

#endif

#pragma once
#include <sqlite3.h>
#include <vector>
#include <string>
#include <functional>


struct ReactionRow {
    int reaction_id;
    int number_of_reactants;
    int number_of_products;
    int reactant_1;
    int reactant_2;
    int product_1;
    int product_2;
    double rate;

    static std::string sql_statement;

    static std::vector<
        std::function<
            void(
                ReactionRow&,
                sqlite3_stmt*,
                int
                )>> getters;

};

std::string ReactionRow::sql_statement =
    "SELECT reaction_id, number_of_reactants, number_of_products, "
    "reactant_1, reactant_2, product_1, product_2, rate FROM reactions;";

std::vector<std::function<
                void(
                    ReactionRow&,
                    sqlite3_stmt*,
                    int)>> ReactionRow::getters = {

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.reaction_id = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.number_of_reactants = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.number_of_products = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.reactant_1 = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.reactant_2 = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.product_1 = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.product_2 = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.rate = sqlite3_column_double(stmt, i);
    }
};

struct TrajectoriesRow {
    int seed;
    int step;
    int reaction_id;
    double time;

    static std::string sql_statement;

    static std::vector<
        std::function<
            int(
                TrajectoriesRow&,
                sqlite3_stmt*,
                int)>> setters;
};

std::string TrajectoriesRow::sql_statement =
    "INSERT INTO trajectories VALUES (?1, ?2, ?3, ?4);";

std::vector<
    std::function<
        int(
            TrajectoriesRow&,
            sqlite3_stmt*,
            int)>> TrajectoriesRow::setters = {

    [](TrajectoriesRow& t, sqlite3_stmt* stmt, int n) {
        return sqlite3_bind_int(stmt, n + 1, t.seed);
    },

    [](TrajectoriesRow& t, sqlite3_stmt* stmt, int n) {
        return sqlite3_bind_int(stmt, n + 1, t.step);
    },

    [](TrajectoriesRow& t, sqlite3_stmt* stmt, int n) {
        return sqlite3_bind_int(stmt, n + 1, t.reaction_id);
    },

    [](TrajectoriesRow& t, sqlite3_stmt* stmt, int n) {
        return sqlite3_bind_double(stmt, n + 1, t.time);
    }
};



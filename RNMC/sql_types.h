#pragma once
#include <sqlite3.h>
#include <vector>
#include <string>
#include <functional>


struct MetadataRow {
    int number_of_species;
    int number_of_reactions;

    static std::string sql_statement;

    static std::vector<
        std::function<
            void(
                MetadataRow&,
                sqlite3_stmt*,
                int
                )>> getters;

};

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


struct FactorsRow {
    double factor_zero;
    double factor_two;
    double factor_duplicate;

    static std::string sql_statement;


    static std::vector<
        std::function<
            void(
                FactorsRow&,
                sqlite3_stmt*,
                int
                )>> getters;


};

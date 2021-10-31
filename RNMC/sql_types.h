#pragma once
#include <sqlite3.h>
#include <vector>
#include <string>
#include <functional>
#include "../core/sql.h"

struct MetadataSql {
    int number_of_species;
    int number_of_reactions;
    static std::string sql_statement;
    static void action(MetadataSql &r, sqlite3_stmt *stmt);
};

struct ReactionSql {
    int reaction_id;
    int number_of_reactants;
    int number_of_products;
    int reactant_1;
    int reactant_2;
    int product_1;
    int product_2;
    double rate;
    static std::string sql_statement;
    static void action(ReactionSql &r, sqlite3_stmt *stmt);
};

struct TrajectoriesSql {
    int seed;
    int step;
    int reaction_id;
    double time;
    static std::string sql_statement;
    static void action(TrajectoriesSql &r, sqlite3_stmt *stmt);
};


struct FactorsSql {
    double factor_zero;
    double factor_two;
    double factor_duplicate;
    static std::string sql_statement;
    static void action(FactorsSql &r, sqlite3_stmt *stmt);
};

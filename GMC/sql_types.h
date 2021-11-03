#pragma once
#include <sqlite3.h>
#include <string>

struct MetadataSql {
    unsigned long int number_of_species;
    unsigned long int number_of_reactions;
    static std::string sql_statement;
    static void action(MetadataSql &r, sqlite3_stmt *stmt);
};

struct ReactionSql {
    unsigned long int reaction_id;
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

struct InitialStateSql {
    int species_id;
    int count;
    static std::string sql_statement;
    static void action(InitialStateSql &r, sqlite3_stmt *stmt);
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


// it is important that these all appear in their own
// compilation unit and not the header. Having them in the
// header will cause different compilation units to have their
// own definitions which will confuse the linker.

std::string MetadataSql::sql_statement =
    "SELECT number_of_species, number_of_reactions FROM metadata;";

void MetadataSql::action(MetadataSql &r, sqlite3_stmt *stmt) {
        r.number_of_species = sqlite3_column_int(stmt, 0);
        r.number_of_reactions = sqlite3_column_int(stmt, 1);
};


std::string ReactionSql::sql_statement =
    "SELECT reaction_id, number_of_reactants, number_of_products, "
    "reactant_1, reactant_2, product_1, product_2, rate FROM reactions;";


void ReactionSql::action(ReactionSql &r, sqlite3_stmt *stmt) {
        r.reaction_id = sqlite3_column_int(stmt, 0);
        r.number_of_reactants = sqlite3_column_int(stmt, 1);
        r.number_of_products = sqlite3_column_int(stmt, 2);
        r.reactant_1 = sqlite3_column_int(stmt, 3);
        r.reactant_2 = sqlite3_column_int(stmt, 4);
        r.product_1 = sqlite3_column_int(stmt, 5);
        r.product_2 = sqlite3_column_int(stmt, 6);
        r.rate = sqlite3_column_double(stmt, 7);
};

std::string InitialStateSql::sql_statement =
    "SELECT species_id, count FROM initial_state;";

void InitialStateSql::action(InitialStateSql &r, sqlite3_stmt *stmt) {
    r.species_id = sqlite3_column_int(stmt, 0);
    r.count = sqlite3_column_int(stmt, 1);
}


std::string TrajectoriesSql::sql_statement =
    "INSERT INTO trajectories VALUES (?1, ?2, ?3, ?4);";

void TrajectoriesSql::action (TrajectoriesSql& t, sqlite3_stmt* stmt) {
    sqlite3_bind_int(stmt, 1, t.seed);
    sqlite3_bind_int(stmt, 2, t.step);
    sqlite3_bind_int(stmt, 3, t.reaction_id);
    sqlite3_bind_double(stmt, 4, t.time);
};

std::string FactorsSql::sql_statement =
    "SELECT factor_zero, factor_two, factor_duplicate FROM factors";


void FactorsSql::action (FactorsSql &r, sqlite3_stmt *stmt) {
    r.factor_zero = sqlite3_column_double(stmt, 0);
    r.factor_two = sqlite3_column_double(stmt, 1);
    r.factor_duplicate = sqlite3_column_double(stmt, 2);
};


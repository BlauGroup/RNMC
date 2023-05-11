#ifndef RNMC_SQL_TYPES_H
#define RNMC_SQL_TYPES_H

#include <sqlite3.h>
#include <string>

struct MetadataSql {
    unsigned long int number_of_species;
    unsigned long int number_of_reactions;
    static std::string sql_statement;
    static void action(MetadataSql &r, sqlite3_stmt *stmt);
};

std::string MetadataSql::sql_statement =
    "SELECT number_of_species, number_of_reactions FROM metadata;";

void MetadataSql::action(MetadataSql &r, sqlite3_stmt *stmt) {
        r.number_of_species = sqlite3_column_int(stmt, 0);
        r.number_of_reactions = sqlite3_column_int(stmt, 1);
};

struct CutoffHistoryElement{
    unsigned long int seed;
    int step;
    double time;
};

/* --------- Read Cutoff SQL ---------*/

struct ReadCutoffSql {
    int seed;
    int step;
    double time;
    static std::string sql_statement;
    static void action(ReadCutoffSql &r, sqlite3_stmt *stmt);
};

std::string ReadCutoffSql::sql_statement =
    "SELECT seed, step, time FROM interrupt_cutoff;";

void ReadCutoffSql::action(ReadCutoffSql &r, sqlite3_stmt *stmt) {
    r.seed = sqlite3_column_int(stmt, 0);
    r.step = sqlite3_column_int(stmt, 1);
    r.time = sqlite3_column_double(stmt, 2);
}

/* --------- WriteCutoff SQL ---------*/

struct WriteCutoffSql {
    int seed;
    int step;
    double time;
    static std::string sql_statement;
    static void action(WriteCutoffSql &r, sqlite3_stmt *stmt);
};

std::string WriteCutoffSql::sql_statement =
    "INSERT INTO interrupt_cutoff VALUES (?1,?2,?3);";

void WriteCutoffSql::action(WriteCutoffSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.step);
    sqlite3_bind_double(stmt, 3, r.time);
}

/* ------------ Factor Sql ------------*/

struct FactorsSql {
    double factor_zero;
    double factor_two;
    double factor_duplicate;
    static std::string sql_statement;
    static void action(FactorsSql &r, sqlite3_stmt *stmt);
};


std::string FactorsSql::sql_statement =
    "SELECT factor_zero, factor_two, factor_duplicate FROM factors";


void FactorsSql::action (FactorsSql &r, sqlite3_stmt *stmt) {
    r.factor_zero = sqlite3_column_double(stmt, 0);
    r.factor_two = sqlite3_column_double(stmt, 1);
    r.factor_duplicate = sqlite3_column_double(stmt, 2);
};

struct InitialStateSql {
    int species_id;
    int count;
    static std::string sql_statement;
    static void action(InitialStateSql &r, sqlite3_stmt *stmt);
};

std::string InitialStateSql::sql_statement =
    "SELECT species_id, count FROM initial_state;";

void InitialStateSql::action(InitialStateSql &r, sqlite3_stmt *stmt) {
    r.species_id = sqlite3_column_int(stmt, 0);
    r.count = sqlite3_column_int(stmt, 1);
}
#endif
#include "sql_types.h"

std::string MetadataSql::sql_statement =
    "SELECT number_of_species, number_of_reactions FROM metadata;";

void MetadataSql::action(MetadataSql &r, sqlite3_stmt *stmt) {
        r.number_of_species = sqlite3_column_int(stmt, 0);
        r.number_of_reactions = sqlite3_column_int(stmt, 1);
}

std::string ReadCutoffSql::sql_statement =
    "SELECT seed, step, time FROM interrupt_cutoff;";

void ReadCutoffSql::action(ReadCutoffSql &r, sqlite3_stmt *stmt) {
    r.seed = sqlite3_column_int(stmt, 0);
    r.step = sqlite3_column_int(stmt, 1);
    r.time = sqlite3_column_double(stmt, 2);
}

std::string WriteCutoffSql::sql_statement =
    "INSERT INTO interrupt_cutoff VALUES (?1,?2,?3);";

void WriteCutoffSql::action(WriteCutoffSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.step);
    sqlite3_bind_double(stmt, 3, r.time);
}

std::string FactorsSql::sql_statement =
    "SELECT factor_zero, factor_two, factor_duplicate FROM factors";


void FactorsSql::action (FactorsSql &r, sqlite3_stmt *stmt) {
    r.factor_zero = sqlite3_column_double(stmt, 0);
    r.factor_two = sqlite3_column_double(stmt, 1);
    r.factor_duplicate = sqlite3_column_double(stmt, 2);
}

std::string InitialStateSql::sql_statement =
    "SELECT species_id, count FROM initial_state;";

void InitialStateSql::action(InitialStateSql &r, sqlite3_stmt *stmt) {
    r.species_id = sqlite3_column_int(stmt, 0);
    r.count = sqlite3_column_int(stmt, 1);
}
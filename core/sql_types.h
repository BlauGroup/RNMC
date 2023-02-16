#pragma once
#include <sqlite3.h>
#include <string>

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
    "SELECT seed, step, time FROM interupt_cutoff;";

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
    "INSERT INTO interupt_cutoff VALUES (?1,?2,?3);";

void WriteCutoffSql::action(WriteCutoffSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.step);
    sqlite3_bind_double(stmt, 3, r.time);
}
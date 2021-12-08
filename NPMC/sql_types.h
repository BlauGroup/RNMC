#pragma once
#include <sqlite3.h>
#include <string>

struct SpecieSql {
    int species_id;
    int degrees_of_freedom;
    static std::string sql_statement;
    static void action(SpecieSql &r, sqlite3_stmt *stmt);
};


struct SiteSql {
    int site_id;
    double x;
    double y;
    double z;
    int species_id;
    static std::string sql_statement;
    static void action(SiteSql &r, sqlite3_stmt *stmt);
};


struct InteractionSql {
    int interaction_id;
    int species_id_1;
    int species_id_2;
    int left_state_1;
    int left_state_2;
    int right_state_1;
    int right_state_2;
    double rate;
    static std::string sql_statement;
    static void action(SiteSql &r, sqlite3_stmt *stmt);
};



#pragma once
#include <sqlite3.h>
#include <string>

struct SpeciesSql {
    int species_id;
    int degrees_of_freedom;
    static std::string sql_statement;
    static void action(SpeciesSql &r, sqlite3_stmt *stmt);
};

std::string SpeciesSql::sql_statement =
    "SELECT species_id, degrees_of_freedom FROM species;";

void SpeciesSql::action(SpeciesSql &r, sqlite3_stmt *stmt) {
    r.species_id = sqlite3_column_int(stmt, 0);
    r.degrees_of_freedom = sqlite3_column_int(stmt, 1);
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

std::string SiteSql::sql_statement =
    "SELECT site_id, x, y, z, species_id FROM sites;";


void SiteSql::action(SiteSql &r, sqlite3_stmt *stmt) {
    r.site_id = sqlite3_column_int(stmt, 0);
    r.x = sqlite3_column_double(stmt, 1);
    r.y = sqlite3_column_double(stmt, 2);
    r.z = sqlite3_column_double(stmt, 3);
    r.species_id = sqlite3_column_int(stmt, 4);
}


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
    static void action(InteractionSql &r, sqlite3_stmt *stmt);
};

std::string InteractionSql::sql_statement =
    "SELECT interaction_id, species_id_1, species_id_2, "
    "left_state_1, left_state_2, right_state_1, right_state_2, "
    "rate FROM intearctions;";

void InteractionSql::action(InteractionSql &r, sqlite3_stmt *stmt) {
    r.interaction_id = sqlite3_column_int(stmt, 0);
    r.species_id_1 = sqlite3_column_int(stmt, 1);
    r.species_id_2 = sqlite3_column_int(stmt, 2);
    r.left_state_1 = sqlite3_column_int(stmt, 3);
    r.left_state_2 = sqlite3_column_int(stmt, 4);
    r.right_state_1 = sqlite3_column_int(stmt, 5);
    r.right_state_2 = sqlite3_column_int(stmt, 6);
    r.rate = sqlite3_column_double(stmt, 7);
}


struct MetadataSql {
    int number_of_species;
    int number_of_sites;
    int number_of_interactions;
    static std::string sql_statement;
    static void action(MetadataSql &r, sqlite3_stmt *stmt);
};

std::string MetadataSql::sql_statement =
    "SELECT number_of_species, number_of_sites, "
    "number_of_interactions FROM metadata;";

void MetadataSql::action(MetadataSql &r, sqlite3_stmt *stmt) {
    r.number_of_species = sqlite3_column_int(stmt, 0);
    r.number_of_sites = sqlite3_column_int(stmt, 1);
    r.number_of_interactions = sqlite3_column_int(stmt, 2);
};

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
    unsigned int species_id;
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
    int number_of_sites;
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
    "SELECT interaction_id, number_of_sites, species_id_1, species_id_2, "
    "left_state_1, left_state_2, right_state_1, right_state_2, "
    "rate FROM interactions;";

void InteractionSql::action(InteractionSql &r, sqlite3_stmt *stmt) {
    r.interaction_id = sqlite3_column_int(stmt, 0);
    r.number_of_sites = sqlite3_column_int(stmt, 1);
    r.species_id_1 = sqlite3_column_int(stmt, 2);
    r.species_id_2 = sqlite3_column_int(stmt, 3);
    r.left_state_1 = sqlite3_column_int(stmt, 4);
    r.left_state_2 = sqlite3_column_int(stmt, 5);
    r.right_state_1 = sqlite3_column_int(stmt, 6);
    r.right_state_2 = sqlite3_column_int(stmt, 7);
    r.rate = sqlite3_column_double(stmt, 8);
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


struct FactorsSql {
    double one_site_interaction_factor;
    double two_site_interaction_factor;
    double interaction_radius_bound;
    std::string distance_factor_type;
    static std::string sql_statement;
    static void action(FactorsSql &r, sqlite3_stmt *stmt);
};

std::string FactorsSql::sql_statement =
    "SELECT one_site_interaction_factor, two_site_interaction_factor, "
    "interaction_radius_bound, distance_factor_type FROM factors;";

void FactorsSql::action(FactorsSql &r, sqlite3_stmt *stmt) {
    r.one_site_interaction_factor = sqlite3_column_double(stmt, 0);
    r.two_site_interaction_factor = sqlite3_column_double(stmt, 1);
    r.interaction_radius_bound = sqlite3_column_double(stmt, 2);

    // generally, casting from raw bytes to char is bad, but downstream,
    // we directly check the resulting string against a finite list of
    // accepted strings, so if garbage gets put in, the program will abort.
    const char *distance_factor_type_raw = (char *) sqlite3_column_text(stmt, 3);
    std::string distance_factor_type ( distance_factor_type_raw );
    r.distance_factor_type = std::move(distance_factor_type);


}

struct InitialStateSql {
    int site_id;
    int degree_of_freedom;
    static std::string sql_statement;
    static void action(InitialStateSql &r, sqlite3_stmt *stmt);
};

std::string InitialStateSql::sql_statement =
    "SELECT site_id, degree_of_freedom FROM initial_state;";


void InitialStateSql::action(InitialStateSql &r, sqlite3_stmt *stmt) {
    r.site_id = sqlite3_column_int(stmt, 0);
    r.degree_of_freedom = sqlite3_column_int(stmt, 1);
}

struct ReadTrajectoriesSql {
    int seed;
    int step;
    double time;
    int site_id_1;
    int site_id_2;
    int interaction_id;
    static std::string sql_statement;
    static void action(ReadTrajectoriesSql &r, sqlite3_stmt *stmt);
};

std::string ReadTrajectoriesSql::sql_statement =
    "SELECT seed, step, time, site_id_1, site_id_2, interaction_id FROM trajectories WHERE seed = ?1;";

void ReadTrajectoriesSql::action(ReadTrajectoriesSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    r.step = sqlite3_column_int(stmt, 2);
    r.time = sqlite3_column_double(stmt, 3);
    r.site_id_1 = sqlite3_column_int(stmt, 4);
    r.site_id_2 = sqlite3_column_int(stmt, 5);
    r.interaction_id = sqlite3_column_int(stmt, 6);
}

struct WriteTrajectoriesSql {
    int seed;
    int step;
    double time;
    int site_id_1;
    int site_id_2;
    int interaction_id;
    static std::string sql_statement;
    static void action(WriteTrajectoriesSql &r, sqlite3_stmt *stmt);
};

std::string WriteTrajectoriesSql::sql_statement =
    "INSERT INTO trajectories VALUES (?1,?2,?3,?4,?5,?6);";

void WriteTrajectoriesSql::action(WriteTrajectoriesSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.step);
    sqlite3_bind_double(stmt, 3, r.time);
    sqlite3_bind_int(stmt, 4, r.site_id_1);
    sqlite3_bind_int(stmt, 5, r.site_id_2);
    sqlite3_bind_int(stmt, 6, r.interaction_id);
}

struct ReadStateSql {
    int seed;
    int site_id;
    int degree_of_freedom;
    static std::string sql_statement;
    static void action(ReadStateSql &r, sqlite3_stmt *stmt);
};

std::string ReadStateSql::sql_statement =
    "SELECT seed, site_id, degree_of_freedom FROM interupt_state;";

void ReadStateSql::action(ReadStateSql &r, sqlite3_stmt *stmt) {
    r.seed = sqlite3_column_int(stmt, 0);
    r.site_id = sqlite3_column_int(stmt, 1);
    r.degree_of_freedom = sqlite3_column_int(stmt, 2);
}

struct WriteStateSql {
    int seed;
    int site_id;
    int degree_of_freedom;
    static std::string sql_statement;
    static void action(WriteStateSql &r, sqlite3_stmt *stmt);
};

std::string WriteStateSql::sql_statement =
    "INSERT INTO interupt_state VALUES (?1,?2,?3);";

void WriteStateSql::action(WriteStateSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.site_id);
    sqlite3_bind_int(stmt, 3, r.degree_of_freedom);
}
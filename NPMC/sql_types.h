#pragma once
#include <sqlite3.h>
#include <string>

/* --------- Species SQL ---------*/
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

/* --------- Site SQL ---------*/

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

/* --------- Interaction SQL ---------*/

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

/* --------- Metadata SQL ---------*/

struct NanoMetadataSql {
    int number_of_species;
    int number_of_sites;
    int number_of_interactions;
    static std::string sql_statement;
    static void action(NanoMetadataSql &r, sqlite3_stmt *stmt);
};

std::string NanoMetadataSql::sql_statement =
    "SELECT number_of_species, number_of_sites, "
    "number_of_interactions FROM metadata;";

void NanoMetadataSql::action(NanoMetadataSql &r, sqlite3_stmt *stmt) {
    r.number_of_species = sqlite3_column_int(stmt, 0);
    r.number_of_sites = sqlite3_column_int(stmt, 1);
    r.number_of_interactions = sqlite3_column_int(stmt, 2);
};

/* --------- Factors SQL ---------*/

struct NanoFactorsSql {
    double one_site_interaction_factor;
    double two_site_interaction_factor;
    double interaction_radius_bound;
    std::string distance_factor_type;
    static std::string sql_statement;
    static void action(NanoFactorsSql &r, sqlite3_stmt *stmt);
};

std::string NanoFactorsSql::sql_statement =
    "SELECT one_site_interaction_factor, two_site_interaction_factor, "
    "interaction_radius_bound, distance_factor_type FROM factors;";

void NanoFactorsSql::action(NanoFactorsSql &r, sqlite3_stmt *stmt) {
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

/* --------- Initial State SQL ---------*/

struct NanoInitialStateSql {
    int site_id;
    int degree_of_freedom;
    static std::string sql_statement;
    static void action(NanoInitialStateSql &r, sqlite3_stmt *stmt);
};

std::string NanoInitialStateSql::sql_statement =
    "SELECT site_id, degree_of_freedom FROM initial_state;";


void NanoInitialStateSql::action(NanoInitialStateSql &r, sqlite3_stmt *stmt) {
    r.site_id = sqlite3_column_int(stmt, 0);
    r.degree_of_freedom = sqlite3_column_int(stmt, 1);
}

/* --------- Read Trajectories SQL ---------*/

struct NanoReadTrajectoriesSql {
    int seed;
    int step;
    double time;
    int site_id_1;
    int site_id_2;
    int interaction_id;
    static std::string sql_statement;
    static void action(NanoReadTrajectoriesSql &r, sqlite3_stmt *stmt);
};

std::string NanoReadTrajectoriesSql::sql_statement =
    "SELECT seed, step, time, site_id_1, site_id_2, interaction_id FROM trajectories;";

void NanoReadTrajectoriesSql::action(NanoReadTrajectoriesSql &r, sqlite3_stmt *stmt) {
    r.seed = sqlite3_column_int(stmt, 0);
    r.step = sqlite3_column_int(stmt, 1);
    r.time = sqlite3_column_double(stmt, 2);
    r.site_id_1 = sqlite3_column_int(stmt, 3);
    r.site_id_2 = sqlite3_column_int(stmt, 4);
    r.interaction_id = sqlite3_column_int(stmt, 5);
}

/* --------- Write Trajectories SQL ---------*/

struct NanoWriteTrajectoriesSql {
    int seed;
    int step;
    double time;
    int site_id_1;
    int site_id_2;
    int interaction_id;
    static std::string sql_statement;
    static void action(NanoWriteTrajectoriesSql &r, sqlite3_stmt *stmt);
};

std::string NanoWriteTrajectoriesSql::sql_statement =
    "INSERT INTO trajectories VALUES (?1,?2,?3,?4,?5,?6);";

void NanoWriteTrajectoriesSql::action(NanoWriteTrajectoriesSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.step);
    sqlite3_bind_double(stmt, 3, r.time);
    sqlite3_bind_int(stmt, 4, r.site_id_1);
    sqlite3_bind_int(stmt, 5, r.site_id_2);
    sqlite3_bind_int(stmt, 6, r.interaction_id);
}

/* --------- Read State SQL ---------*/

struct NanoReadStateSql {
    int seed;
    int site_id;
    int degree_of_freedom;
    static std::string sql_statement;
    static void action(NanoReadStateSql &r, sqlite3_stmt *stmt);
};

std::string NanoReadStateSql::sql_statement =
    "SELECT seed, site_id, degree_of_freedom FROM interrupt_state;";

void NanoReadStateSql::action(NanoReadStateSql &r, sqlite3_stmt *stmt) {
    r.seed = sqlite3_column_int(stmt, 0);
    r.site_id = sqlite3_column_int(stmt, 1);
    r.degree_of_freedom = sqlite3_column_int(stmt, 2);
}

/* --------- Write State SQL ---------*/
struct NanoWriteStateSql {
    int seed;
    int site_id;
    int degree_of_freedom;
    static std::string sql_statement;
    static void action(NanoWriteStateSql &r, sqlite3_stmt *stmt);
};

std::string NanoWriteStateSql::sql_statement =
    "INSERT INTO interrupt_state VALUES (?1,?2,?3);";

void NanoWriteStateSql::action(NanoWriteStateSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.site_id);
    sqlite3_bind_int(stmt, 3, r.degree_of_freedom);
}

/* --------- State and Trajectory History Elements ---------*/
// Each Element will be stored in a History Packet which will be stored 
// in a history queue to be dumped to SQL database in batches

struct NanoStateHistoryElement{
    unsigned long int seed; //seed
    int site_id; 
    int degree_of_freedom; //energy level the site is at
};

struct NanoTrajectoryHistoryElement {

    unsigned long int seed; // seed
    NanoReaction reaction; // reaction which fired
    double time;  // time after reaction has occoured.
    int step;
};
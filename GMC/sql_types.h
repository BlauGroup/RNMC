#pragma once
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

/* ------------ Write Trajectory ------------*/

struct ReactionNetworkWriteTrajectoriesSql {
    int seed;
    int step;
    int reaction_id;
    double time;
    static std::string sql_statement;
    static void action(ReactionNetworkWriteTrajectoriesSql &r, sqlite3_stmt *stmt);
};

std::string ReactionNetworkWriteTrajectoriesSql::sql_statement =
    "INSERT INTO trajectories VALUES (?1, ?2, ?3, ?4);";

void ReactionNetworkWriteTrajectoriesSql::action (ReactionNetworkWriteTrajectoriesSql& t, sqlite3_stmt* stmt) {
    sqlite3_bind_int(stmt, 1, t.seed);
    sqlite3_bind_int(stmt, 2, t.step);
    sqlite3_bind_int(stmt, 3, t.time);
    sqlite3_bind_double(stmt, 4, t.reaction_id);
};

/* ------------ Read Trajectory ------------*/

struct ReactionNetworkReadTrajectoriesSql {
    int seed;
    int step;
    double time;
    int reaction_id;
    static std::string sql_statement;
    static void action(ReactionNetworkReadTrajectoriesSql &r, sqlite3_stmt *stmt);
};

std::string ReactionNetworkReadTrajectoriesSql::sql_statement =
    "SELECT seed, step, time, reaction_id FROM trajectories;";

void ReactionNetworkReadTrajectoriesSql::action(ReactionNetworkReadTrajectoriesSql &r, sqlite3_stmt *stmt) {
    r.seed = sqlite3_column_int(stmt, 0);
    r.step = sqlite3_column_int(stmt, 1);
    r.time = sqlite3_column_double(stmt, 2);
    r.reaction_id = sqlite3_column_int(stmt, 3);
    
}
/* ------------ Read state ------------*/

struct ReactionNetworkReadStateSql {
    int seed;
    int species_id;
    int count;
    static std::string sql_statement;
    static void action(ReactionNetworkReadStateSql &r, sqlite3_stmt *stmt);
};

std::string ReactionNetworkReadStateSql::sql_statement =
    "SELECT seed, species_id, count FROM interrupt_state;";

void ReactionNetworkReadStateSql::action(ReactionNetworkReadStateSql &r, sqlite3_stmt *stmt) {
    r.seed = sqlite3_column_int(stmt, 0);
    r.species_id = sqlite3_column_int(stmt, 1);
    r.count = sqlite3_column_int(stmt, 2);
    
}

/* ------------ Write state ------------*/

struct ReactionNetworkWriteStateSql {
    int seed;
    int species_id;
    int count;
    static std::string sql_statement;
    static void action(ReactionNetworkWriteStateSql &r, sqlite3_stmt *stmt);
};

std::string ReactionNetworkWriteStateSql::sql_statement =
    "INSERT INTO interrupt_state VALUES (?1,?2,?3);";

void ReactionNetworkWriteStateSql::action(ReactionNetworkWriteStateSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.species_id);
    sqlite3_bind_int(stmt, 3, r.count);
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

/* --------- State and Trajectory History Elements ---------*/
// Each Element will be stored in a History Packet which will be stored 
// in a history queue to be dumped to SQL database in batches

struct ReactionNetworkStateHistoryElement{
    unsigned long int seed; //seed
    int species_id;
    int count;
};

struct ReactionNetworkTrajectoryHistoryElement {

    unsigned long int seed; // seed
    int reaction_id; // reaction which fired
    double time;  // time after reaction has occoured.
    int step;
};
#pragma once
#include <sqlite3.h>
#include <string>

// TODO: Dealing with strings
// Should these v char's be strings or char*'s?

struct LGMCReactionSql {
    unsigned long int reaction_id;
    
    int number_of_reactants;
    int number_of_products;
    
    int reactant_1;
    int reactant_2;
    int product_1;
    int product_2;
    
    const unsigned char* phase_reactant_1;
    const unsigned char* phase_reactant_2;
    const unsigned char* phase_product_1;
    const unsigned char* phase_product_2;
    
    double dG;
    double prefactor;
    double rate;

    double electron_tunneling_coefficient;
    double reorganization_energy;
    double charge_transfer_coefficient;

    const unsigned char* type;
    
    static std::string sql_statement;
    static void action(LGMCReactionSql &r, sqlite3_stmt *stmt);
};

std::string LGMCReactionSql::sql_statement =
    "SELECT reaction_id, number_of_reactants, number_of_products, "
    "reactant_1, reactant_2, product_1, product_2, phase_reactant_1, phase_reactant_2, "
    "phase_product_1, phase_product_2, dG, prefactor, rate, electron_tunneling_coefficient, "
    "reorganization_energy, charge_transfer_coefficient, type FROM reactions;";


void LGMCReactionSql::action(LGMCReactionSql &r, sqlite3_stmt *stmt) {
        r.reaction_id = sqlite3_column_int(stmt, 0);
        r.number_of_reactants = sqlite3_column_int(stmt, 1);
        r.number_of_products = sqlite3_column_int(stmt, 2);
        r.reactant_1 = sqlite3_column_int(stmt, 3);
        r.reactant_2 = sqlite3_column_int(stmt, 4);
        r.product_1 = sqlite3_column_int(stmt, 5);
        r.product_2 = sqlite3_column_int(stmt, 6);
        r.phase_reactant_1 = sqlite3_column_text(stmt, 7);
        r.phase_reactant_2 = sqlite3_column_text(stmt, 8);
        r.phase_product_1 = sqlite3_column_text(stmt, 9);
        r.phase_product_2 = sqlite3_column_text(stmt, 10);
        r.dG = sqlite3_column_double(stmt, 11);
        r.prefactor = sqlite3_column_double(stmt, 12);
        r.rate = sqlite3_column_double(stmt, 13);
        r.electron_tunneling_coefficient = sqlite3_column_double(stmt, 14);
        r.reorganization_energy = sqlite3_column_double(stmt, 15);
        r.charge_transfer_coefficient = sqlite3_column_double(stmt, 16);
        r.type = sqlite3_column_text(stmt, 17);
};


//TODO test and deal with text

/* --------- Read Trajectories SQL ---------*/

struct LatticeReadTrajectoriesSql {
    int seed;
    int step;
    double time;
    int reaction_id;
    int site_1;
    int site_2;
    static std::string sql_statement;
    static void action(LatticeReadTrajectoriesSql &r, sqlite3_stmt *stmt);
};

std::string LatticeReadTrajectoriesSql::sql_statement =
    "SELECT seed, step, time, reaction_id, site_1, site_2 FROM trajectories;";

void LatticeReadTrajectoriesSql::action (LatticeReadTrajectoriesSql& r, sqlite3_stmt* stmt) {
    r.seed = sqlite3_column_int(stmt, 0);
    r.step = sqlite3_column_int(stmt, 1);
    r.time = sqlite3_column_double(stmt, 2);
    r.reaction_id = sqlite3_column_int(stmt, 3);
    r.site_1 = sqlite3_column_int(stmt, 4);
    r.site_2 = sqlite3_column_int(stmt, 5);
    
};

/* --------- Write Trajectories SQL ---------*/

struct LatticeWriteTrajectoriesSql {
    int seed;
    int step;
    double time;
    int reaction_id;
    int site_1;
    int site_2;
    static std::string sql_statement;
    static void action(LatticeWriteTrajectoriesSql &r, sqlite3_stmt *stmt);
};

std::string LatticeWriteTrajectoriesSql::sql_statement =
    "INSERT INTO trajectories VALUES (?1, ?2, ?3, ?4, ?5, ?6);";

void LatticeWriteTrajectoriesSql::action (LatticeWriteTrajectoriesSql& t, sqlite3_stmt* stmt) {
    sqlite3_bind_int(stmt, 1, t.seed);
    sqlite3_bind_int(stmt, 2, t.step);
    sqlite3_bind_double(stmt, 3, t.time);
    sqlite3_bind_int(stmt, 4, t.reaction_id);
    sqlite3_bind_int(stmt, 5, t.site_1);
    sqlite3_bind_int(stmt, 6, t.site_2);
    
};

/* --------- Read State SQL ---------*/

struct LatticeReadStateSql {
    int seed;
    int site_id;
    int species_id;
    int quantity;
    static std::string sql_statement;
    static void action(LatticeReadStateSql &r, sqlite3_stmt *stmt);
};

std::string LatticeReadStateSql::sql_statement =
    "SELECT seed, site_1, site_2, reaction_id FROM interupt_state;";

void LatticeReadStateSql::action(LatticeReadStateSql &r, sqlite3_stmt *stmt) {
    r.seed = sqlite3_column_int(stmt, 0);
    r.site_id = sqlite3_column_int(stmt, 1);
    r.species_id = sqlite3_column_int(stmt, 2);
    r.quantity = sqlite3_column_int(stmt, 3);
}

/* --------- Write State SQL ---------*/
struct LatticeWriteStateSql {
    int seed;
    int site_id;
    int species_id;
    int quantity;
    static std::string sql_statement;
    static void action(LatticeWriteStateSql &r, sqlite3_stmt *stmt);
};

std::string LatticeWriteStateSql::sql_statement =
    "INSERT INTO interupt_state VALUES (?1,?2,?3, ?4);";

void LatticeWriteStateSql::action(LatticeWriteStateSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.site_id);
    sqlite3_bind_int(stmt, 3, r.species_id);
    sqlite3_bind_int(stmt, 4, r.quantity);
}

/* --------- State and Trajectory History Elements ---------*/

struct LatticeStateHistoryElement{
    unsigned long int seed; //seed
    int site_id;
    int species_id;
    int quantity;
};

struct LatticeTrajectoryHistoryElement {

    unsigned long int seed; // seed
    int step;
    double time;  // time after reaction has occoured.
    int reaction_id; // reaction which fired
    int site_1;
    int site_2;

};

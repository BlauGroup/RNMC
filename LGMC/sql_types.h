#ifndef RNMC_LGMC_SQL_TYPES_H
#define RNMC_LGMC_SQL_TYPES_H

#include <sqlite3.h>
#include <string>

class LGMCReactionSql {
public:
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

/* --------- Read Trajectories SQL ---------*/

class LatticeReadTrajectoriesSql {
public: 
    int seed;
    int step;
    double time;
    int reaction_id;
    int site_1_mapping;
    int site_2_mapping;
    static std::string sql_statement;
    static void action(LatticeReadTrajectoriesSql &r, sqlite3_stmt *stmt);
};

/* --------- Write Trajectories SQL ---------*/

class LatticeWriteTrajectoriesSql {
public:
    int seed;
    int step;
    double time;
    int reaction_id;
    int site_1_mapping;
    int site_2_mapping;
    static std::string sql_statement;
    static void action(LatticeWriteTrajectoriesSql &r, sqlite3_stmt *stmt);
};

/* --------- Read State SQL ---------*/

class LatticeReadStateSql {
public:
    int seed;
    int species_id;
    int quantity;
    int site_mapping;
    int edge;
    static std::string sql_statement;
    static void action(LatticeReadStateSql &r, sqlite3_stmt *stmt);
};

/* --------- Write State SQL ---------*/

class LatticeWriteStateSql {
public:
    int seed;
    int species_id;
    int quantity;
    int site_mapping;
    int edge;
    static std::string sql_statement;
    static void action(LatticeWriteStateSql &r, sqlite3_stmt *stmt);
};

/* --------- Read Cutoff SQL ---------*/

class LatticeReadCutoffSql {
public:
    int seed;
    int step;
    double time;
    int maxk;
    static std::string sql_statement;
    static void action(LatticeReadCutoffSql &r, sqlite3_stmt *stmt);
};

/* --------- WriteCutoff SQL ---------*/

class LatticeWriteCutoffSql {
public:
    int seed;
    int step;
    double time;
    int maxk;
    static std::string sql_statement;
    static void action(LatticeWriteCutoffSql &r, sqlite3_stmt *stmt);
};

/* --------- State, Trajectory, and Cutoff History Elements ---------*/

struct LatticeStateHistoryElement{
    unsigned long int seed; 
    int species_id;
    int quantity;
    int site_mapping;
    int edge;
};

struct LatticeTrajectoryHistoryElement {

    unsigned long int seed; 
    int step;
    double time;  
    int reaction_id; 
    int site_1_mapping;
    int site_2_mapping;

};

struct LatticeCutoffHistoryElement{
    unsigned long int seed;
    int step;
    double time;
    int maxk;
};


#endif
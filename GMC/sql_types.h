#ifndef RNMC_GMC_SQL_TYPES_H
#define RNMC_GMC_SQL_TYPES_H

#include <sqlite3.h>
#include <string>

class ReactionSql {
public:
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

class ReactionNetworkWriteTrajectoriesSql {
public:
    int seed;
    int step;
    int reaction_id;
    double time;
    static std::string sql_statement;
    static void action(ReactionNetworkWriteTrajectoriesSql &r, sqlite3_stmt *stmt);
};

class ReactionNetworkReadTrajectoriesSql {
public:
    int seed;
    int step;
    double time;
    int reaction_id;
    static std::string sql_statement;
    static void action(ReactionNetworkReadTrajectoriesSql &r, sqlite3_stmt *stmt);
};

class ReactionNetworkReadStateSql {
public:
    int seed;
    int species_id;
    int count;
    static std::string sql_statement;
    static void action(ReactionNetworkReadStateSql &r, sqlite3_stmt *stmt);
};

class ReactionNetworkWriteStateSql {
public:
    int seed;
    int species_id;
    int count;
    static std::string sql_statement;
    static void action(ReactionNetworkWriteStateSql &r, sqlite3_stmt *stmt);
};

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

#endif
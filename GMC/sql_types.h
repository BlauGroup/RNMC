/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#ifndef RNMC_GMC_SQL_TYPES_H
#define RNMC_GMC_SQL_TYPES_H

#include <sqlite3.h>
#include <string>

class ReactionSql
{
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

class EnergyReactionSql
{
public:
    unsigned long int reaction_id;
    int number_of_reactants;
    int number_of_products;
    int reactant_1;
    int reactant_2;
    int product_1;
    int product_2;
    double rate;
    double dG;
    static std::string sql_statement;
    static void action(EnergyReactionSql &r, sqlite3_stmt *stmt);
};

/* --------- I/O Trajectories SQL ---------*/

class ReactionNetworkWriteTrajectoriesSql
{
public:
    int seed;
    int step;
    int reaction_id;
    double time;
    static std::string sql_statement;
    static void action(ReactionNetworkWriteTrajectoriesSql &r, sqlite3_stmt *stmt);
};

class ReactionNetworkReadTrajectoriesSql
{
public:
    int seed;
    int step;
    double time;
    int reaction_id;
    static std::string sql_statement;
    static void action(ReactionNetworkReadTrajectoriesSql &r, sqlite3_stmt *stmt);
};

/* --------- I/O State SQL ---------*/

class ReactionNetworkReadStateSql
{
public:
    int seed;
    int species_id;
    int count;
    static std::string sql_statement;
    static void action(ReactionNetworkReadStateSql &r, sqlite3_stmt *stmt);
};

class ReactionNetworkWriteStateSql
{
public:
    int seed;
    int species_id;
    int count;
    static std::string sql_statement;
    static void action(ReactionNetworkWriteStateSql &r, sqlite3_stmt *stmt);
};

/* --------- I/O Cutoff SQL ---------*/

class EnergyNetworkReadCutoffSql
{
public:
    int seed;
    int step;
    double time;
    double energy_budget;
    static std::string sql_statement;
    static void action(EnergyNetworkReadCutoffSql &r, sqlite3_stmt *stmt);
};

class EnergyNetworkWriteCutoffSql
{
public:
    int seed;
    int step;
    double time;
    double energy_budget;
    static std::string sql_statement;
    static void action(EnergyNetworkWriteCutoffSql &r, sqlite3_stmt *stmt);
};

/* --------- State, Trajectory, and Cutoff History Elements ---------*/

struct ReactionNetworkStateHistoryElement
{
    unsigned long int seed; 
    int species_id;
    int count;
};

struct ReactionNetworkTrajectoryHistoryElement
{
    unsigned long int seed; 
    int reaction_id;        // reaction which fired
    double time;            // time after reaction has occoured.
    int step;
};

struct EnergyNetworkCutoffHistoryElement
{
    unsigned long int seed;
    int step;
    double time;
    double energy_budget;
};

#endif
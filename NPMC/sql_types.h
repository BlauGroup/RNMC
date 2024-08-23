/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#ifndef RNMC_NPMC_SQL_TYPES_H
#define RNMC_NPMC_SQL_TYPES_H

#include <sqlite3.h>
#include <string>

#include "NPMC_types.h"

/* --------------------------- Species SQL --------------------------- */

class SpeciesSql
{
public:
    int species_id;
    int degrees_of_freedom;
    static std::string sql_statement;
    static void action(SpeciesSql &r, sqlite3_stmt *stmt);
};

/* --------------------------- Site SQL ------------------------------ */

class SiteSql
{
public:
    int site_id;
    double x;
    double y;
    double z;
    unsigned int species_id;
    static std::string sql_statement;
    static void action(SiteSql &r, sqlite3_stmt *stmt);
};

/* ----------------------- Interaction SQL --------------------------- */

class InteractionSql
{
public:
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

/* ---------------------------- Metadata SQL ------------------------- */

class NanoMetadataSql
{
public:
    int number_of_species;
    int number_of_sites;
    int number_of_interactions;
    static std::string sql_statement;
    static void action(NanoMetadataSql &r, sqlite3_stmt *stmt);
};

/* ----------------------------- Factors SQL ------------------------- */

class NanoFactorsSql
{
public:
    double one_site_interaction_factor;
    double two_site_interaction_factor;
    double interaction_radius_bound;
    std::string distance_factor_type;
    static std::string sql_statement;
    static void action(NanoFactorsSql &r, sqlite3_stmt *stmt);
};

/* ------------------------- Initial state SQL ----------------------- */


class NanoInitialStateSql
{
public:
    int site_id;
    int degree_of_freedom;
    static std::string sql_statement;
    static void action(NanoInitialStateSql &r, sqlite3_stmt *stmt);
};

/* ------------------------ I/O Trajectories SQL --------------------- */

class NanoReadTrajectoriesSql
{
public:
    int seed;
    int step;
    double time;
    int site_id_1;
    int site_id_2;
    int interaction_id;
    static std::string sql_statement;
    static void action(NanoReadTrajectoriesSql &r, sqlite3_stmt *stmt);
};

class NanoWriteTrajectoriesSql
{
public:
    int seed;
    int step;
    double time;
    int site_id_1;
    int site_id_2;
    int interaction_id;
    static std::string sql_statement;
    static void action(NanoWriteTrajectoriesSql &r, sqlite3_stmt *stmt);
};

/* ----------------------------- I/O state SQL ----------------------- */

class NanoReadStateSql
{
public:
    int seed;
    int site_id;
    int degree_of_freedom;
    static std::string sql_statement;
    static void action(NanoReadStateSql &r, sqlite3_stmt *stmt);
};

class NanoWriteStateSql
{
public:
    int seed;
    int site_id;
    int degree_of_freedom;
    static std::string sql_statement;
    static void action(NanoWriteStateSql &r, sqlite3_stmt *stmt);
};

/* ----------------- State and Trajectory History Elements ----------- */

struct NanoStateHistoryElement
{
    unsigned long int seed; // seed
    int site_id;
    int degree_of_freedom; // energy level the site is at
};

struct NanoTrajectoryHistoryElement
{
    unsigned long int seed; // seed
    NanoReaction reaction;  // reaction which fired
    double time;            // time after reaction has occoured.
    int step;
};

#endif
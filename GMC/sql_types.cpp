/* ----------------------------------------------------------------------
RNMC - Reaction Network Monte Carlo
https://blaugroup.github.io/RNMC/

See the README file in the top-level RNMC directory.
---------------------------------------------------------------------- */

#include "sql_types.h"

/* -------------------------------- Reaction --------------------------------*/

std::string ReactionSql::sql_statement =
    "SELECT reaction_id, number_of_reactants, number_of_products, "
    "reactant_1, reactant_2, product_1, product_2, rate FROM reactions;";

void ReactionSql::action(ReactionSql &r, sqlite3_stmt *stmt)
{
    r.reaction_id = sqlite3_column_int(stmt, 0);
    r.number_of_reactants = sqlite3_column_int(stmt, 1);
    r.number_of_products = sqlite3_column_int(stmt, 2);
    r.reactant_1 = sqlite3_column_int(stmt, 3);
    r.reactant_2 = sqlite3_column_int(stmt, 4);
    r.product_1 = sqlite3_column_int(stmt, 5);
    r.product_2 = sqlite3_column_int(stmt, 6);
    r.rate = sqlite3_column_double(stmt, 7);
};

/* ----------------------------- Energy Reaction -----------------------------*/


std::string EnergyReactionSql::sql_statement =
    "SELECT reaction_id, number_of_reactants, number_of_products, "
    "reactant_1, reactant_2, product_1, product_2, rate, dG FROM reactions;";

void EnergyReactionSql::action(EnergyReactionSql &r, sqlite3_stmt *stmt)
{
    r.reaction_id = sqlite3_column_int(stmt, 0);
    r.number_of_reactants = sqlite3_column_int(stmt, 1);
    r.number_of_products = sqlite3_column_int(stmt, 2);
    r.reactant_1 = sqlite3_column_int(stmt, 3);
    r.reactant_2 = sqlite3_column_int(stmt, 4);
    r.product_1 = sqlite3_column_int(stmt, 5);
    r.product_2 = sqlite3_column_int(stmt, 6);
    r.rate = sqlite3_column_double(stmt, 7);
    r.dG = sqlite3_column_double(stmt, 8);
};

/* ----------------------------- Write Trajectory -----------------------------*/

std::string ReactionNetworkWriteTrajectoriesSql::sql_statement =
    "INSERT INTO trajectories VALUES (?1, ?2, ?3, ?4);";

void ReactionNetworkWriteTrajectoriesSql::action(ReactionNetworkWriteTrajectoriesSql &t, sqlite3_stmt *stmt)
{
    sqlite3_bind_int(stmt, 1, t.seed);
    sqlite3_bind_int(stmt, 2, t.step);
    sqlite3_bind_int(stmt, 3, t.reaction_id);
    sqlite3_bind_double(stmt, 4, t.time);
};

/* ----------------------------- Read Trajectory -----------------------------*/

std::string ReactionNetworkReadTrajectoriesSql::sql_statement =
    "SELECT seed, step, reaction_id, time FROM trajectories;";

void ReactionNetworkReadTrajectoriesSql::action(ReactionNetworkReadTrajectoriesSql &r, sqlite3_stmt *stmt)
{
    r.seed = sqlite3_column_int(stmt, 0);
    r.step = sqlite3_column_int(stmt, 1);
    r.reaction_id = sqlite3_column_int(stmt, 2);
    r.time = sqlite3_column_double(stmt, 3);
}

/* --------------------------------- Read state ------------------------------*/

std::string ReactionNetworkReadStateSql::sql_statement =
    "SELECT seed, species_id, count FROM interrupt_state;";

void ReactionNetworkReadStateSql::action(ReactionNetworkReadStateSql &r, sqlite3_stmt *stmt)
{
    r.seed = sqlite3_column_int(stmt, 0);
    r.species_id = sqlite3_column_int(stmt, 1);
    r.count = sqlite3_column_int(stmt, 2);
}

/* ------------------------------ Write state --------------------------------*/

std::string ReactionNetworkWriteStateSql::sql_statement =
    "INSERT INTO interrupt_state VALUES (?1,?2,?3);";

void ReactionNetworkWriteStateSql::action(ReactionNetworkWriteStateSql &r, sqlite3_stmt *stmt)
{
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.species_id);
    sqlite3_bind_int(stmt, 3, r.count);
}

/* ---------------- State and Trajectory History Elements ----------------------*/
// Each Element will be stored in a History Packet which will be stored
// in a history queue to be dumped to SQL database in batches

std::string EnergyNetworkReadCutoffSql::sql_statement =
    "SELECT seed, step, time, energy_budget FROM interrupt_cutoff;";

void EnergyNetworkReadCutoffSql::action(EnergyNetworkReadCutoffSql &r, sqlite3_stmt *stmt)
{
    r.seed = sqlite3_column_int(stmt, 0);
    r.step = sqlite3_column_int(stmt, 1);
    r.time = sqlite3_column_double(stmt, 2);
    r.energy_budget = sqlite3_column_double(stmt, 3);
}

std::string EnergyNetworkWriteCutoffSql::sql_statement =
    "INSERT INTO interrupt_cutoff VALUES (?1,?2,?3, ?4);";

void EnergyNetworkWriteCutoffSql::action(EnergyNetworkWriteCutoffSql &r, sqlite3_stmt *stmt)
{
    sqlite3_bind_int(stmt, 1, r.seed);
    sqlite3_bind_int(stmt, 2, r.step);
    sqlite3_bind_double(stmt, 3, r.time);
    sqlite3_bind_double(stmt, 4, r.energy_budget);
}
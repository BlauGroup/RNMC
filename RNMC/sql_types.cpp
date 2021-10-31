#include "sql_types.h"

// it is important that these all appear in their own
// compilation unit and not the header. Having them in the
// header will cause different compilation units to have their
// own definitions which will confuse the linker.

std::string MetadataRow::sql_statement =
    "SELECT number_of_species, number_of_reactions FROM metadata;";

void MetadataRow::action(MetadataRow &r, sqlite3_stmt *stmt) {
        r.number_of_species = sqlite3_column_int(stmt, 0);
        r.number_of_reactions = sqlite3_column_int(stmt, 1);
};


std::string ReactionRow::sql_statement =
    "SELECT reaction_id, number_of_reactants, number_of_products, "
    "reactant_1, reactant_2, product_1, product_2, rate FROM reactions;";


void ReactionRow::action(ReactionRow &r, sqlite3_stmt *stmt) {
        r.reaction_id = sqlite3_column_int(stmt, 0);
        r.number_of_reactants = sqlite3_column_int(stmt, 1);
        r.number_of_products = sqlite3_column_int(stmt, 2);
        r.reactant_1 = sqlite3_column_int(stmt, 3);
        r.reactant_2 = sqlite3_column_int(stmt, 4);
        r.product_1 = sqlite3_column_int(stmt, 5);
        r.product_2 = sqlite3_column_int(stmt, 6);
        r.rate = sqlite3_column_double(stmt, 7);
};

std::string TrajectoriesRow::sql_statement =
    "INSERT INTO trajectories VALUES (?1, ?2, ?3, ?4);";

void TrajectoriesRow::action (TrajectoriesRow& t, sqlite3_stmt* stmt) {
    sqlite3_bind_int(stmt, 1, t.seed);
    sqlite3_bind_int(stmt, 2, t.step);
    sqlite3_bind_int(stmt, 3, t.reaction_id);
    sqlite3_bind_double(stmt, 4, t.time);
};

std::string FactorsRow::sql_statement =
    "SELECT factor_zero, factor_two, factor_duplicate FROM factors";


void FactorsRow::action (FactorsRow &r, sqlite3_stmt *stmt) {
    r.factor_zero = sqlite3_column_double(stmt, 1);
    r.factor_two = sqlite3_column_double(stmt, 2);
    r.factor_duplicate = sqlite3_column_double(stmt, 3);

};


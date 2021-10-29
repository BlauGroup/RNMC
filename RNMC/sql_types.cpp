#include "sql_types.h"

// it is important that these all appear in their own
// compilation unit and not the header. Having them in the
// header will cause different compilation units to have their
// own definitions which will confuse the linker.

std::string MetadataRow::sql_statement =
    "SELECT number_of_species, number_of_reactions FROM metadata;";

std::vector<
    std::function<
        void(
            MetadataRow&,
            sqlite3_stmt*,
            int
            )>> MetadataRow::getters = {

    [](MetadataRow &r, sqlite3_stmt *stmt, int i) {
        r.number_of_species = sqlite3_column_int(stmt, i);
    },

    [](MetadataRow &r, sqlite3_stmt *stmt, int i) {
        r.number_of_reactions = sqlite3_column_int(stmt, i);
    },
};

std::string ReactionRow::sql_statement =
    "SELECT reaction_id, number_of_reactants, number_of_products, "
    "reactant_1, reactant_2, product_1, product_2, rate FROM reactions;";

std::vector<std::function<
                void(
                    ReactionRow&,
                    sqlite3_stmt*,
                    int)>> ReactionRow::getters = {

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.reaction_id = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.number_of_reactants = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.number_of_products = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.reactant_1 = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.reactant_2 = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.product_1 = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.product_2 = sqlite3_column_int(stmt, i);
    },

    [](ReactionRow &r, sqlite3_stmt *stmt, int i) {
        r.rate = sqlite3_column_double(stmt, i);
    }
};

std::string TrajectoriesRow::sql_statement =
    "INSERT INTO trajectories VALUES (?1, ?2, ?3, ?4);";

std::vector<
    std::function<
        int(
            TrajectoriesRow&,
            sqlite3_stmt*,
            int)>> TrajectoriesRow::setters = {

    [](TrajectoriesRow& t, sqlite3_stmt* stmt, int n) {
        return sqlite3_bind_int(stmt, n + 1, t.seed);
    },

    [](TrajectoriesRow& t, sqlite3_stmt* stmt, int n) {
        return sqlite3_bind_int(stmt, n + 1, t.step);
    },

    [](TrajectoriesRow& t, sqlite3_stmt* stmt, int n) {
        return sqlite3_bind_int(stmt, n + 1, t.reaction_id);
    },

    [](TrajectoriesRow& t, sqlite3_stmt* stmt, int n) {
        return sqlite3_bind_double(stmt, n + 1, t.time);
    }
};

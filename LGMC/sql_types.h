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

struct LGMCTrajectoriesSql {
    int seed;
    int step;
    int reaction_id;
    int site_1;
    int site_2;
    double time;
    static std::string sql_statement;
    static void action(LGMCTrajectoriesSql &r, sqlite3_stmt *stmt);
};

std::string LGMCTrajectoriesSql::sql_statement =
    "INSERT INTO trajectories VALUES (?1, ?2, ?3, ?4, ?5, ?6);";

void LGMCTrajectoriesSql::action (LGMCTrajectoriesSql& t, sqlite3_stmt* stmt) {
    sqlite3_bind_int(stmt, 1, t.seed);
    sqlite3_bind_int(stmt, 2, t.step);
    sqlite3_bind_int(stmt, 3, t.reaction_id);
    sqlite3_bind_int(stmt, 4, t.site_1);
    sqlite3_bind_int(stmt, 5, t.site_2);
    sqlite3_bind_double(stmt, 6, t.time);
};
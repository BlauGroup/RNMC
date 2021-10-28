#include <sqlite3.h>
#include <string>
#include <utility>
#include <vector>
#include <functional>


class SqlConnection {
public:
    sqlite3 *connection;
    std::string database_file_path;

    // method for executing standalone sql statements.
    // for reading and writing data, use SqlReader or SqlWriter classes.
    void exec(std::string sql_statement) {
        sqlite3_exec(
            connection,
            sql_statement.c_str(),
            nullptr,
            nullptr,
            nullptr);
    };

    SqlConnection(std::string database_file_path) :
        database_file_path{database_file_path} {
            sqlite3_open(database_file_path.c_str(), &connection);
    };

    ~SqlConnection() {
        sqlite3_close(connection);
    };

    // no copy constructor because we don't have access to internal state
    // of a sql connection.
    SqlConnection(SqlConnection &other) = delete;

    // move constructor
    SqlConnection(SqlConnection &&other) :
        connection{std::exchange(other.connection, nullptr)},
        database_file_path{std::move(other.database_file_path)} {};

    // no copy assignment because we don't have access to internal state
    // of a sql connection.
    SqlConnection &operator=(SqlConnection &other) = delete;

    // move assignment
    SqlConnection &operator=(SqlConnection &&other) {
        std::swap(connection, other.connection);
        std::swap(database_file_path, other.database_file_path);
        return *this;
    };
};

template<typename T>
// T needs sql_statement, setters, nothing and empty attributes.
class SqlReader {
private:
    sqlite3_stmt *stmt;
    SqlConnection &sql_connection_ref;
public:
    // important that readers and writers have a reference to
    // the database connection. This is because it is an error
    // to close a connection before all statements have been
    // finalized.

    T next() {
        int rc = sqlite3_step(stmt);

        if (rc != SQLITE_OK) {
            return T::empty;
        }

        T result;
        result.nothing = false;

        for (int i = 0; i < T::setters.size(); i++) {
            T::setters[i](*result, sql_connection_ref.connection, i);
        }

        return result;
    };

    SqlReader(SqlConnection &sql_connection_ref) :

        // TODO: deal with error handling
        // not a huge issue as every sql statement is hard coded
        // into a Row type.
        sql_connection_ref{sql_connection_ref} {
            int rc = sqlite3_prepare_v2(
                sql_connection_ref.connection,
                T::sql_statement.c_str(),
                -1,
                &stmt,
                nullptr
                );
    };

    ~SqlReader() {
        sqlite3_finalize(stmt);
    }

    // no copy constructor
    // don't have access to sqlite internals so can't copy sqlite statements
    SqlReader(SqlReader &other) = delete;

    // move constructor
    SqlReader(SqlReader &&other) :
        sqlite3_stmt{std::exchange(other.stmt, nullptr)},
        sql_connection_ref{other.sql_connection_ref} {
    };

    // no copy assigment
    // don't have access to sqlite internals so can't copy sqlite statements
    SqlReader &operator=(SqlReader &other) = delete;

    // move assignment
    SqlReader &operator=(SqlReader &&other) {
        std::swap(stmt, other.stmt);
        std::swap(sql_connection_ref, other.sql_connection_ref);
        return *this;
    };
};

// Row structs correspond to rows in a sqlite database.
// The setters attribute is a static vector of functions which
// we can call to set the corresponding attributes in the Row struct
// according to the column numbers from the sql statement.

struct ReactionRow {
    bool nothing;
    int reaction_id;
    int number_of_reactants;
    int number_of_products;
    int reactant_1;
    int reactant_2;
    int product_1;
    int product_2;
    double rate;

    static std::string sql_statement;

    static std::vector<
        std::function<
            void(
                ReactionRow&,
                sqlite3_stmt*,
                int
                )>> setters;

    static ReactionRow empty;
};

ReactionRow empty = ReactionRow {
    .nothing = true,
    .reaction_id = -1,
    .number_of_reactants = -1,
    .number_of_products = -1,
    .reactant_1 = -1,
    .reactant_2 = -1,
    .product_1 = -1,
    .product_2 = -1,
    .rate = -1.0};

std::string sql_statement =
    "SELECT reaction_id, number_of_reactants, number_of_products, "
    "reactant_1, reactant_2, product_1, product_2, rate FROM reactions;";

std::vector<std::function<
                void(
                    ReactionRow&,
                    sqlite3_stmt*,
                    int)>> ReactionRow::setters = {

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

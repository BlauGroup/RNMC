#pragma once
#include <sqlite3.h>
#include <string>
#include <utility>
#include <vector>
#include <functional>
#include <optional>
#include <iostream>

// IO API used in all the simulators.  light wrapper around sqlite.
// key concept is row structs Row structs correspond to rows in a
// sqlite database.  The getters attribute is a static vector of
// functions which we can call to get the corresponding attributes in
// the Row struct The setters attribute is a static vector of
// functions which we can call to set the corresponding attributes in
// the Row struct according to the column numbers from the sql
// statement.  recall, in SQL, ?n variables are 1 indexed, which is
// why all these lambdas have n + 1 in them. Here are some examples.

// struct ExampleSelectRow {
//     int a;
//     double b;

//     static std::string sql_statement;

//     static std::vector<
//         std::function<
//             void(
//                 ExampleSelectRow&,
//                 sqlite3_stmt*,
//                 int
//                 )>> getters;
// };

// std::string ExampleSelectRow::sql_statement = "SELECT a, b FROM foo;";

// std::vector<std::function<
//                 void(
//                     ExampleSelectRow&,
//                     sqlite3_stmt*,
//                     int)>> ExampleSelectRow::getters = {

//     [](ExampleSelectRow &r, sqlite3_stmt *stmt, int i) {
//         r.a = sqlite3_column_int(stmt, i);
//     },

//     [](ExampleSelectRow &r, sqlite3_stmt *stmt, int i) {
//         r.b = sqlite3_column_double(stmt, i);
//     }
// };

// struct ExampleInsertRow {
//     int a;
//     double b;

//     static std::string sql_statement;

//     static std::vector<
//         std::function<
//             int(
//                 ExampleInsertRow&,
//                 sqlite3_stmt*,
//                 int)>> setters;
// };

// std::string ExampleInsertRow::sql_statement =
//     "INSERT INTO trajectories VALUES (?1, ?2);";

// std::vector<
//     std::function<
//         int(
//             ExampleInsertRow&,
//             sqlite3_stmt*,
//             int)>> ExampleInsertRow::setters = {

//     [](ExampleInsertRow& r, sqlite3_stmt* stmt, int n) {
//         return sqlite3_bind_int(stmt, n + 1, r.a);
//     },

//     [](ExampleInsertRow& r, sqlite3_stmt* stmt, int n) {
//         return sqlite3_bind_double(stmt, n + 1, r.b);
//     }
// };



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
        database_file_path (database_file_path) {
            int rc = sqlite3_open_v2(
                database_file_path.c_str(),
                &connection,
                SQLITE_OPEN_READWRITE,
                nullptr
                );

            if (rc != SQLITE_OK) {
                std::cerr << "sqlite: "
                          << sqlite3_errmsg(connection)
                          << '\n';
                std::abort();
            }
    };

    ~SqlConnection() {
        sqlite3_close(connection);
    };

    // no copy constructor because we don't have access to internal state
    // of a sql connection.
    SqlConnection(SqlConnection &other) = delete;

    // move constructor
    SqlConnection(SqlConnection &&other) :
        connection (std::exchange(other.connection, nullptr)),
        database_file_path (std::move(other.database_file_path)) {};

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
// T needs sql_statement and getters attributes.
class SqlReader {
private:
    sqlite3_stmt *stmt;
    SqlConnection &sql_connection;
    bool done;
public:
    // important that readers and writers have a reference to
    // the database connection. This is because it is an error
    // to close a connection before all statements have been
    // finalized.

    std::optional<T> next() {
        if (done) return std::optional<T> ();
        else {

            int rc = sqlite3_step(stmt);

            if (rc == SQLITE_DONE) {
                done = true;
                return std::optional<T> ();
            }

            T result;

            for (int i = 0; i < T::getters.size(); i++) {
                T::getters[i](std::ref(result), stmt, i);
            }

            return std::optional<T> (result);
        };
    };

    // TODO: write a reset method so we can reloop without needing to create
    // a new object.

    SqlReader(SqlConnection &sql_connection) :

        sql_connection (sql_connection),
        done (false) {
            int rc = sqlite3_prepare_v2(
                sql_connection.connection,
                T::sql_statement.c_str(),
                -1,
                &stmt,
                nullptr
                );

            if (rc != SQLITE_OK) {
                std::cerr << "sqlite: "
                          << sqlite3_errmsg(sql_connection.connection)
                          << '\n';

                std::abort();
            }
    };

    ~SqlReader() {
        sqlite3_finalize(stmt);
    }

    // no copy constructor
    // don't have access to sqlite internals so can't copy sqlite statements
    SqlReader(SqlReader &other) = delete;

    // move constructor
    SqlReader(SqlReader &&other) :
        sqlite3_stmt (std::exchange(other.stmt, nullptr)),
        sql_connection (other.sql_connection),
        done (other.done)
        {
    };

    // no copy assigment
    // don't have access to sqlite internals so can't copy sqlite statements
    SqlReader &operator=(SqlReader &other) = delete;

    // move assignment
    SqlReader &operator=(SqlReader &&other) {
        std::swap(stmt, other.stmt);
        std::swap(sql_connection, other.sql_connection);
        done = other.done;
        return *this;
    };
};

template<typename T>
// T needs sql_statement and setters attributes.
class SqlWriter {
private:
    sqlite3_stmt *stmt;
    SqlConnection &sql_connection;

public:
    SqlWriter(SqlConnection &sql_connection) :
        sql_connection (sql_connection)
        {
            int rc = sqlite3_prepare_v2(
                sql_connection.connection,
                T::sql_statement.c_str(),
                -1,
                &stmt,
                nullptr
                );

            if (rc != SQLITE_OK) {
                std::cerr << "sqlite: "
                          << sqlite3_errmsg(sql_connection.connection)
                          << '\n';

                std::abort();
            }
        };

    ~SqlWriter() {
        sqlite3_finalize(stmt);
    }

    void insert(T row) {
        sqlite3_reset(stmt);
        for (int i = 0; i < T::setters.size(); i++) {
            T::setters[i](row, stmt, i);
        }

        // TODO: error handling
        int rc = sqlite3_step(stmt);
    };

    // no copy constructor
    // don't have access to sqlite internals so can't copy sqlite statements
    SqlWriter(SqlWriter &other) = delete;

    // move constructor
    SqlWriter(SqlWriter &&other) :
        sqlite3_stmt (std::exchange(other.stmt, nullptr)),
        sql_connection(other.sql_connection) {};

    // no copy assigment
    // don't have access to sqlite internals so can't copy sqlite statements
    SqlWriter &operator=(SqlWriter &other) = delete;

    // move assignment
    SqlWriter &operator=(SqlWriter &&other) {
        std::swap(stmt, other.stmt);
        std::swap(sql_connection, other.sql_connection);
        return *this;
    };

};

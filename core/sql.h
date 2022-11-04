#pragma once
#include <sqlite3.h>
#include <string>
#include <vector>
#include <optional>
#include <iostream>
#include <iomanip>
#include <cstring>

// DESIGN
// the sqlite C API has two kinds of resources: connections and
// statements. We store connections in the class SqlConnection. We store
// statements in the class SqlStatement which is parameterized over a
// SqlType. Two example SqlTypes are shown below. They consist of a
// bunch of attributes, a static sql_statement and a static
// action. The action is the C code required to link the attributes to
// the sql statement.

auto time_stamp() {
    auto time = std::time(nullptr);
    return std::put_time(std::localtime(&time), "[%T] ");
};

struct ExampleSelectSql {
    int foo;
    double bar;
    static std::string sql_statement;
    static void action(ExampleSelectSql &r, sqlite3_stmt *stmt);
};


struct ExampleInsertSql {
    int foo;
    double bar;
    static std::string sql_statement;
    static void action(ExampleInsertSql &r, sqlite3_stmt *stmt);
};

std::string ExampleSelectSql::sql_statement = "SELECT foo, bar FROM table;";
void ExampleSelectSql::action(ExampleSelectSql &r, sqlite3_stmt *stmt) {
    r.foo = sqlite3_column_int(stmt, 0);
    r.bar = sqlite3_column_double(stmt, 1);
}

std::string ExampleInsertSql::sql_statement = "INSERT INTO table VALUES (?1, ?2);";
void ExampleInsertSql::action(ExampleInsertSql &r, sqlite3_stmt *stmt) {
    sqlite3_bind_int(stmt, 1, r.foo);
    sqlite3_bind_double(stmt, 2, r.bar);
}

class SqlConnection {
public:
    sqlite3 *connection;
    std::string database_file_path;

    // method for executing standalone sql statements.
    // for reading and writing data, use SqlReader or SqlWriter classes.
    void exec(std::string sql_statement) {
        int rc;
        char *error_message = 0;
        
        rc = sqlite3_exec(
            connection,
            sql_statement.c_str(),
            nullptr,
            nullptr,
            &error_message);
        if( rc != SQLITE_OK ){
            char error[200];
            strcpy(error,"SQL error: ");
            strcat(error,error_message);
            std::cerr << error << std::endl;
        }
    };

    void close() {
        sqlite3_close(connection);
    }
    SqlConnection(std::string database_file_path, int sql_flags) :
        database_file_path (database_file_path) {
            int rc = sqlite3_open_v2(
                database_file_path.c_str(),
                &connection,
                sql_flags,
                nullptr
                );

            if (rc != SQLITE_OK) {
                std::cerr << time_stamp()
                          << "sqlite: "
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



template <typename T>
class SqlStatement {
private:

    // important that SqlStatement objects have a reference to
    // the database connection. This is because it is an error
    // to close a connection before all statements have been
    // finalized.
    sqlite3_stmt *stmt;
    SqlConnection &sql_connection;

public:
    void action(T &r) { T::action(r, stmt);};
    void reset() { sqlite3_reset(stmt); };
    int step() { return sqlite3_step(stmt); };


    SqlStatement(SqlConnection &sql_connection) :
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
                std::cerr << time_stamp()
                          << "sqlite: "
                          << sqlite3_errmsg(sql_connection.connection)
                          << '\n';

                std::abort();
            }
        };

    ~SqlStatement() {
        sqlite3_finalize(stmt);
    }

    // no copy constructor
    // don't have access to sqlite internals so can't copy sqlite statements
    SqlStatement(SqlStatement &other) = delete;

    // move constructor
    SqlStatement(SqlStatement &&other) :
        sqlite3_stmt (std::exchange(other.stmt, nullptr)),
        sql_connection(other.sql_connection) {};

    // no copy assigment
    // don't have access to sqlite internals so can't copy sqlite statements
    SqlStatement &operator=(SqlStatement &other) = delete;

    // move assignment
    SqlStatement &operator=(SqlStatement &&other) {
        std::swap(stmt, other.stmt);
        std::swap(sql_connection, other.sql_connection);
        return *this;
    };
};



template <typename T>
class SqlReader {
private:
    bool done;
    SqlStatement<T> &statement;

public:

    SqlReader(SqlStatement<T> &statement) :
        done (false),
        statement (statement)
        {};


    // TODO: write a reset method so we can reloop without needing to
    // create a new object
    std::optional<T> next() {
        if (done) {
            return std::optional<T> ();

        } else {

            int rc = statement.step();

            if (rc == SQLITE_DONE) {
                done = true;
                return std::optional<T> ();
            }

            T result;
            statement.action(result);

            return std::optional<T> (result);
        }
    };
};



template <typename T>
class SqlWriter {
private:
    SqlStatement<T> &statement;


public:
    SqlWriter(SqlStatement<T> &statement) :
        statement (statement) {};

    void insert(T row) {
        statement.reset();
        statement.action(row);
        statement.step();

        // TODO: error handling
    };


};

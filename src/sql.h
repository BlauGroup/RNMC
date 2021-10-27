#include <sqlite3.h>
#include <string>
#include <utility>


class SqlConnection {
private:
    sqlite3 *connection;
public:
    std::string database_file_path;

    // method for executing standalone sql statements.
    // for reading and writing data, use SqlReader or SqlWriter classes.
    void exec(std::string statement) {
        sqlite3_exec(
            connection,
            statement.c_str(),
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

// sqlite3_column template. generic case is not defined.
// Fallback definitions follow for specific types.

template<typename T>
T sqlite3_column(sqlite3_stmt *stmt, int col);

template<>
int sqlite3_column<int>(sqlite3_stmt *stmt, int col) {
    return sqlite3_column_int(stmt, col);
}

template<>
double sqlite3_column<double>(sqlite3_stmt *stmt, int col) {
    return sqlite3_column_double(stmt, col);
}



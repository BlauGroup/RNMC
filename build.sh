mkdir -p build

flags="-fno-rtti -fno-exceptions -fno-unwind-tables -std=c++17 -g $(gsl-config --cflags) $(gsl-config --libs) -lsqlite3 -lpthread"

$CC $flags ./core/*.cpp -o ./build/core_tests

mkdir -p build

flags="-O0 -fno-rtti -fno-exceptions -fno-unwind-tables -std=c++17 -g $(gsl-config --cflags) $(gsl-config --libs) -lsqlite3 -lpthread"

$CC $flags ./core/*.cpp -o ./build/core_test

$CC $flags ./RNMC/*.cpp -o ./build/RNMC

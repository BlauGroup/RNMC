mkdir -p build

flags="-fno-rtti -fno-exceptions -std=c++17 -Wall -Wextra -g $(gsl-config --cflags) $(gsl-config --libs) -lsqlite3 -lpthread"

# $CC $flags ./core/test.cpp -o ./build/core_test

$CC $flags ./GMC/GMC.cpp -o ./build/GMC

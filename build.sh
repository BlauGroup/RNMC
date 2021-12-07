mkdir -p build

flags="-fno-rtti -fno-exceptions -std=c++17 -Wall -Wextra -g $(gsl-config --cflags) $(gsl-config --libs) -lsqlite3 -lpthread"

$CC $flags ./core/test.cpp -o ./build/test_core
$CC $flags ./GMC/GMC.cpp -o ./build/GMC
$CC $flags ./NPMC/NPMC.cpp -o ./build/NPMC

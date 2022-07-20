mkdir -p build

flags="-fno-rtti -fno-exceptions -std=c++20 -Wall -Wextra -g $(gsl-config --cflags) $(gsl-config --libs) -lsqlite3 -lpthread"

#echo "building test_core"
#$CC $flags ./core/test.cpp -o ./build/test_core
#echo "building GMC"
#$CC $flags ./GMC/GMC.cpp -o ./build/GMC
#echo "building NPMC"
#$CC $flags ./NPMC/NPMC.cpp -o ./build/NPMC

$CC $flags ./LGMC/LGMC.cpp -o ./build/LGMC

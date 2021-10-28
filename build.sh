$CC -fno-rtti -fno-exceptions -fno-unwind-tables -std=c++17 -g ./src/*.cpp -o NPMC $(gsl-config --cflags) $(gsl-config --libs) -lsqlite3 -lpthread

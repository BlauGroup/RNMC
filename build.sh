$CC -fno-rtti -fno-exceptions -fno-unwind-tables -g ./src/*.cpp -o NPMC $(gsl-config --cflags) $(gsl-config --libs) -lsqlite3 -lpthread

$CC ./src/*.cpp -O0 -g -o NPMC $(gsl-config --cflags) $(gsl-config --libs) -lsqlite3 -lpthread

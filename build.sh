$CC ./src/*.cpp -g -o NPMC $(gsl-config --cflags) $(gsl-config --libs) -lsqlite3 -lpthread

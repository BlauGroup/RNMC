# Installation

## Dependencies
RNMC depends on [GSL](https://www.gnu.org/software/gsl/) for pseudo random number generation and [sqlite](https://www.sqlite.org/index.html) for the database interfaces.

## Makefile

Inside RNMC is a main makefile which can be used for GMC, NPMC, LGMC, or testing. Alternatively there are individual makefiles inside each folder. 

For GMC, NPMC and LGMC:

```$ make GMC/NPMC/LGMC```

```$ make clean```


On a machine with system versions of GSL and sqlite, the executables can be built like this:
```
./build.sh
```
The executables are put in the `build` directory. Note that the build script uses the `gsl-config` utility to find headers and libraries for GSL. If you are on a cluster and sqlite is not present, it can be built as follows:

```
cd $HOME
wget https://www.sqlite.org/2021/sqlite-amalgamation-3360000.zip
unzip sqlite-amalgamation-3360000.zip
cd sqlite-amalgamation-3360000
gcc -o libsqlite3.so -shared -fPIC sqlite3.c -lpthread -ldl
```

in which case, the simulators can be built like this:

```
export CPATH=$HOME/sqlite-amalgamation-3360000:$CPATH
export LIBRARY_PATH=$HOME/sqlite-amalgamation-3360000:$LIBRARY_PATH
./build.sh
```
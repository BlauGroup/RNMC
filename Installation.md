# Installation

## Dependencies
RNMC depends on [GSL](https://www.gnu.org/software/gsl/) for pseudo random number generation and [sqlite](https://www.sqlite.org/index.html) for the database interfaces.

## Makefile
On a machine with system versions of GSL and sqlite, the executables can be built with a makefile. There are makefiles inside the `GMC`, `NPMC`, or `LGMC` folders.

To make an executable for `GMC`, `NPMC` or `LGMC` first enter that folder to use the makefile. To create an executable, use make and then the name of the module. 

For example, to use `GMC`:

```
$ cd GMC
```

```
$ make GMC
```

For further help on the makefile and to view other commands:

```
$ make help
```

Note that the makefile uses the `gsl-config` utility to find headers and libraries for GSL. If you are on a cluster and sqlite is not present, it can be built as follows:

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
make GMC
```

If you need to build [GSL](https://www.gnu.org/software/gsl/) from source: 

```
wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
mkdir gsl
mv gsl-latest.tar.gz gsl
tar -xvf gsl-latest.tar.gz
cd gsl-2.7.1
./configure --prefix=(desired location of gsl)
make
make install
echo $PKG_CONFIG_PATH 
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:(location of gsl-2.7.1)
```

Note that if you build from source use `pkg-config gsl` instead of `gsl-config` inside each makefile

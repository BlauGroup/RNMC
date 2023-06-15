# How to use this Makefile... make help

EXECUTABLE  = NPMC

core_DIR = ../core

# designate which compiler to use
CXX         = g++

# list of sources used
SOURCES     = $(wildcard *.cpp) $(core_DIR)/sql_types.cpp

# list of objects used
OBJECTS     = $(SOURCES:%.cpp=%.o) 

# name of the perf data file, only used by the clean target
PERF_FILE = perf.data*

# Default Flags (we prefer -std=c++17 but Mac/Xcode/Clang doesn't support)
# WARNING: Adding flags like _GLIBCXX_DEBUG or -fsanitize
# may prevent your project from working properly!
CXXFLAGS = -fno-rtti -fno-exceptions -std=c++20 -Wall -Wextra $(shell gsl-config --cflags) $(shell gsl-config --libs) -lsqlite3 -lpthread

# make debug - will compile sources with $(CXXFLAGS) and the -g3 flag
#              also defines DEBUG, so "#ifdef DEBUG /*...*/ #endif" works
debug: CXXFLAGS += -g3 -DDEBUG
debug:
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(EXECUTABLE)_debug 

# make profile - will compile "all" with $(CXXFLAGS) and the -g3 and -O3 flags
profile: CXXFLAGS += -g3 -O3
profile:
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(EXECUTABLE)_profile 
# make gprof - will compile "all" with $(CXXFLAGS) and the -pg (for gprof)
gprof: CXXFLAGS += -pg
gprof:
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(EXECUTABLE)_profile 

# Build all executables
all: debug profile gprof

$(EXECUTABLE): $(OBJECTS)
ifeq ($(EXECUTABLE), executable)
	@echo Edit EXECUTABLE variable in Makefile.
	@echo Using default a.out.
	$(CXX) $(CXXFLAGS) $(OBJECTS) 
else
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)
endif

# rule for creating objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $(dir $@)$(notdir $@)

# make clean - remove .o files, executables, tarball
clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(EXECUTABLE)_debug $(EXECUTABLE)_profile 
	rm -Rf *.dSYM

define MAKEFILE_HELP
Makefile Help
* For more information:
    $$ man make

* General usage
	1. To make an executable for NPMC $$ make NPMC
	2. You can $$ make debug, profile, gprof, or clean
	3. To do all at once $$ make all 

endef
export MAKEFILE_HELP

help:
	@echo "$$MAKEFILE_HELP"

# these targets do not create any files
.PHONY: all debug profile gprof clean 

# disable built-in rules
.SUFFIXES:
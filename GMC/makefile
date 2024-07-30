# How to use this Makefile... make help

EXECUTABLE  = GMC

core_DIR = ../core

# store executable in ./build directory
BUILD_DIR = ../build

# designate which compiler to use
CXX         = g++

# list of sources used
SOURCES     = $(wildcard *.cpp) $(core_DIR)/sql_types.cpp

# list of objects used
OBJECTS     = $(SOURCES:%.cpp=%.o) 

# name of the perf data file, only used by the clean target
PERF_FILE = perf.data*

# flags
CXXFLAGS = -fno-rtti -fno-exceptions -std=c++17 -Wall -Wextra -lsqlite3 -lpthread \
			$(shell gsl-config --cflags) $(shell gsl-config --libs) 

# make debug - will compile sources with $(CXXFLAGS) and the -g3 flag
#              also defines DEBUG, so "#ifdef DEBUG /*...*/ #endif" works
debug: CXXFLAGS += -g3 -DDEBUG
debug:
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(BUILD_DIR)$(EXECUTABLE)_debug 

# make profile - will compile "all" with $(CXXFLAGS) and the -g3 and -O3 flags
profile: CXXFLAGS += -g3 -O3
profile:
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(BUILD_DIR)/$(EXECUTABLE)_profile 
# make gprof - will compile "all" with $(CXXFLAGS) and the -pg (for gprof)
gprof: CXXFLAGS += -pg
gprof:
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(BUILD_DIR)/$(EXECUTABLE)_profile 

# Build all executables
all: debug profile gprof

$(EXECUTABLE): 
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(BUILD_DIR)/$(EXECUTABLE)

# rule for creating objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $(BUILD_DIR)/$(dir $@)$(notdir $@)

# make clean - remove .o files, executables, tarball
clean:
	rm -f $(BUILD_DIR)/$(OBJECTS) $(BUILD_DIR)/$(EXECUTABLE) $(BUILD_DIR)/$(EXECUTABLE)_debug $(BUILD_DIR)/$(EXECUTABLE)_profile 
	rm -Rf *.dSYM

define MAKEFILE_HELP
Makefile Help
* For more information:
    $$ man make

* General usage
	1. To make an executable for GMC $$ make GMC
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

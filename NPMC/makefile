# How to use this Makefile...
###################
###################
##               ##
##  $ make help  ##
##               ##
###################
###################
EXECUTABLE  = NPMC

# designate which compiler to use
CXX         = g++

# list of test drivers (with main()) for development
TESTSOURCES = $(wildcard *test.cpp)
# names of test executables
TESTS       = $(TESTSOURCES:%.cpp=%)

# list of sources used
SOURCES     = $(wildcard *.cpp)
SOURCES     := $(filter-out $(TESTSOURCES), $(SOURCES))

# list of objects used
OBJECTS     = $(SOURCES:%.cpp=%.o)

# name of the perf data file, only used by the clean target
PERF_FILE = perf.data*

# Default Flags (we prefer -std=c++17 but Mac/Xcode/Clang doesn't support)
# WARNING: Adding flags like _GLIBCXX_DEBUG or -fsanitize
# may prevent your project from working properly!
CXXFLAGS = -fno-rtti -fno-exceptions -std=c++20 -Wall -Wextra -g $(shell gsl-config --cflags) $(shell gsl-config --libs) -lsqlite3 -lpthread

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
	$(CXX) $(CXXFLAGS) -c $*.cpp

# make clean - remove .o files, executables, tarball
clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(EXECUTABLE)_debug $(EXECUTABLE)_profile \
      $(TESTS) 
	rm -Rf *.dSYM

define MAKEFILE_HELP
Makefile Help
* This Makefile uses advanced techniques, for more information:
    $$ man make

* General usage
    1. Follow directions at each "TODO" in this file.
       a. Set EXECUTABLE equal to the name from the project specification.
    2. Build, test ... repeat as necessary.

endef
export MAKEFILE_HELP

help:
	@echo "$$MAKEFILE_HELP"

# these targets do not create any files
.PHONY: all debug profile gprof clean 

# disable built-in rules
.SUFFIXES:
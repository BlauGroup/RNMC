# How to use this Makefile...
###################
###################
##               ##
##  $ make help  ##
##               ##
###################
###################
EXECUTABLE  = main

# The following line looks for a project's main() in files named project*.cpp,
# executable.cpp (substituted from EXECUTABLE above), or main.cpp
#PROJECTFILE = $(or $(wildcard project*.cpp $(EXECUTABLE).cpp), main.cpp)
PROJECTFILE = main.cpp

# designate which compiler to use
CXX         = g++

# list of test drivers (with main()) for development
TESTSOURCES = $(wildcard test*.cpp)
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
CXXFLAGS = -fno-rtti -fno-exceptions -std=c++17 -Wall -Wextra -g $(shell gsl-config --cflags) $(shell gsl-config --libs) -lsqlite3 -lpthread

# make release - will compile "all" with $(CXXFLAGS) and the -O3 flag
#                also defines NDEBUG so that asserts will not check
release: CXXFLAGS += -O3 -DNDEBUG
release: $(EXECUTABLE)

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
all:
	+$(MAKE) -C core
	+$(MAKE) -C GMC
	+$(MAKE) -C LGMC
	+$(MAKE) -C NPMC

$(EXECUTABLE): $(OBJECTS)
ifeq ($(EXECUTABLE), executable)
	@echo Edit EXECUTABLE variable in Makefile.
	@echo Using default a.out.
	$(CXX) $(CXXFLAGS) $(OBJECTS)
else
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXECUTABLE)
endif

# Automatically generate any build rules for test*.cpp files
define make_tests
	ifeq ($$(PROJECTFILE),)
		@echo Edit PROJECTFILE variable to .cpp file with main\(\)
		@exit 1
	endif
	SRCS = $$(filter-out $$(PROJECTFILE), $$(SOURCES))
	OBJS = $$(SRCS:%.cpp=%.o)
	HDRS = $$(wildcard *.h *.hpp)
	$(1): CXXFLAGS += -g3 -DDEBUG
	$(1): $$(OBJS) $$(HDRS) $(1).cpp
	$$(CXX) $$(CXXFLAGS) $$(OBJS) $(1).cpp -o $(1)
endef
$(foreach test, $(TESTS), $(eval $(call make_tests, $(test))))

alltests: $(TESTS)

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

#######################
# TODO (begin) #
#######################
# individual dependencies for objects
# Examples:
# "Add a header file dependency"
# project2.o: project2.cpp project2.h
#
# "Add multiple headers and a separate class"
# HEADERS = some.h special.h header.h files.h
# myclass.o: myclass.cpp myclass.h $(HEADERS)
# project5.o: project5.cpp myclass.o $(HEADERS)
#
# SOME EXAMPLES
#
#test_thing: test_thing.cpp class.o functions.o
#class.o: class.cpp class.h
#functions.o: functions.cpp functions.h
#project0.o: project0.cpp class.h functions.h
#
# THE COMPILER CAN GENERATE DEPENDENCIES FROM SOURCE CODE
#
# % g++ -std=c++1z -MM *.cpp
#
# ADD YOUR OWN DEPENDENCIES HERE

######################
# TODO (end) #
######################

# these targets do not create any files
.PHONY: all release debug profile gprof clean 

# disable built-in rules
.SUFFIXES:
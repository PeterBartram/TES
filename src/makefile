CPC = g++
CC = gcc

# C++ compiler flags
CPFLAGS = -Wall -O3 -std=c++14 -static-libstdc++ -fext-numeric-literals -DGIT_REVISION="\"NotWorking\""

# C compiler flags
CFLAGS = -Wall -O3  -DGIT_REVISION="\"NotWorking\"" 

# Module directory
MODULE_DIR = .

# Test directory location
TEST_DIR = ../../TestProjects

INC_DIR = -I$(MODULE_DIR)  \


DEPS = Simulation.c  dhem.c radau.c main.cpp terminate.c radau_step.c UniversalVars.c
DEPS_LONG_DOUBLE = Simulation.c  dhem.c radau.c main.cpp terminate.c radau_step.c universal_variables_ld.c
DEPS_KEPLER = Simulation.c UniversalVars.c


all:
	make build_deps
	make build_problem

long:
	make build_deps_long_double
	make build_long_double

	
build_problem:
	$(CPC) $(CPFLAGS) -o output $(DEPS)  -lpthread -lm 

build_deps:
	$(CC) -c $(DEPS) $(CFLAGS) $(INC_DIR)
	
build_long_double:
	$(CPC) $(CPFLAGS) -DUSE_LONG_DOUBLES -o output $(DEPS_LONG_DOUBLE)  -lpthread -lm

build_deps_long_double:
	$(CC) -c -DUSE_LONG_DOUBLES $(DEPS_LONG_DOUBLE) $(CFLAGS) $(INC_DIR) 	

# Only build the Kepler solver components - used for testing in isolation.
build_kepler:
	$(CC) -c $(DEPS_KEPLER) $(CFLAGS) $(INC_DIR)

clean:
	rm output ||:
	rm *.o ||:

rebuild_deps:
	make clean
	make build_deps

problem:
	make build_deps
	make build_problem

run:
	./output

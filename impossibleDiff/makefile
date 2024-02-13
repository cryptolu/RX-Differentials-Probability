
CXX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++17 -fopenmp -Wno-deprecated-copy #-g -fsanitize=address -fno-omit-frame-pointer

IFLAGS = -I $(GUROBI_HOME)/include/
LFLAGS = -L $(GUROBI_HOME)/lib/ -lgurobi_g++5.2 -lgurobi95

%.o: %.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CXXFLAGS) -c $< -o $@

main : main.o common.o MILP_common.o speck.o customCallback.o lea.o alzette.o xtea.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o main main.o common.o MILP_common.o speck.o customCallback.o lea.o alzette.o xtea.o $(LFLAGS)

allTests : MILP_tests.o common.o MILP_common.o speck.o customCallback.o lea.o alzette.o xtea.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o allTests MILP_tests.o common.o MILP_common.o speck.o customCallback.o lea.o alzette.o xtea.o $(LFLAGS)

testShiftRXDiff : testShiftRXDiff.o common.o
	$(CXX) $(CXXFLAGS) -o testShiftRXDiff testShiftRXDiff.o common.o

clean :
	rm -rf *.o

CXX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++17 -fopenmp -Wno-deprecated-copy #-g -fsanitize=address -fno-omit-frame-pointer

IFLAGS = -I $(GUROBI_HOME)/include/
LFLAGS = -L $(GUROBI_HOME)/lib/ -lgurobi_g++5.2 -lgurobi100

%.o: %.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CXXFLAGS) -c $< -o $@

main : main.o common.o MILP_common.o alzette.o speck.o customCallback.o salsa.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o main main.o common.o MILP_common.o alzette.o speck.o customCallback.o salsa.o $(LFLAGS)

clean :
	rm -rf *.o
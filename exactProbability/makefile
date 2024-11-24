CXX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++17 -fopenmp -Wno-deprecated-copy

IFLAGS = -I $(GUROBI_HOME)/include/

LFLAGS = -L $(GUROBI_HOME)/lib/
LFLAGS_GUROBI = -L $(GUROBI_HOME)/lib/ -lgurobi_g++8.5 -lgurobi110

%.o: %.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CXXFLAGS) -c $< -o $@

common/rxdp.o: common/rxdp.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CXXFLAGS) -c common/rxdp.cpp -o common/rxdp.o $(LFLAGS)

common/customCallback.o: common/customCallback.cpp common/rxdp.o common/common.o
	$(CXX) $(IFLAGS) $(LFLAGS) $(CXXFLAGS) -c common/customCallback.cpp -o common/customCallback.o $(LFLAGS_GUROBI)

verifyTheory: verifyTheory.o common/common.o common/rxdp.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o verifyTheory verifyTheory.o common/rxdp.o common/common.o $(LFLAGS)

runAlzette: runAlzette.o common/common.o common/rxdp.o ciphers/alzette.o common/MILP_common.o common/customCallback.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o runAlzette runAlzette.o common/rxdp.o common/common.o common/MILP_common.o ciphers/alzette.o common/customCallback.o $(LFLAGS_GUROBI)

verifyAlzette: verifyAlzette.o common/common.o common/MILP_common.o common/customCallback.o common/rxdp.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o verifyAlzette verifyAlzette.o common/common.o common/rxdp.o common/MILP_common.o common/customCallback.o $(LFLAGS_GUROBI)

plots:
	bash -xc 'cd logs; for f in *.txt; do python ../process_log_bounds_evolution.py "$f" "../figs/$f.pdf"; done'
	bash -xc 'cd logs; for f in *.txt; do python ../process_log_bounds_evolution.py "$f" "../figs/$f.png"; done'


runA1:
	# c0..c7 k=1..16
	bash -c 'parallel -j 2 ./runAlzette ::: {0..7} ::: {1..16} ::: 0 ::: 3600 ::: 14400' # 1 hour NEQopt
	bash -c 'parallel -j 2 ./runAlzette ::: 0 ::: {1..16} ::: 1 ::: 3600 ::: 14400' # 4 hours FullOpt for c0
	bash -c 'parallel -j 2 ./runAlzette ::: {1..7} ::: {1..16} ::: 1 ::: 3600 ::: 7200' # 2 hours FullOpt for c1..c7


verifyA1:
	MORE_CONSTANTS=0 time ./verifyAlzette trails/alzette_constc2b3293d_rotation2_fullopt_time7200.txt 42 123 456\
		| tee logs/verify_alzette_constc2b3293d_rotation2_fullopt_time7200.txt
	MORE_CONSTANTS=0 time ./verifyAlzette trails/alzette_constc2b3293d_rotation1_fullopt_time7200.txt 42 123 456\
		| tee logs/verify_alzette_constc2b3293d_rotation1_fullopt_time7200.txt

	MORE_CONSTANTS=1 time ./verifyAlzette trails/alzette_constc2b3293d_rotation2_fullopt_time7200.txt 42 123 456\
		| tee logs/verify_alzette_constc2b3293d_rotation2_fullopt_time7200.txt
	MORE_CONSTANTS=1 time ./verifyAlzette trails/alzette_constc2b3293d_rotation1_fullopt_time7200.txt 42 123 456\
		| tee logs/verify_alzette_constc2b3293d_rotation1_fullopt_time7200.txt

verifyA2:
	MORE_CONSTANTS=0 time ./verifyAlzette trails/alzette_constb7e15162_rotation1_fullopt_time14400.txt 42 123 456\
		| tee logs/verify_alzette_constb7e15162_rotation1_fullopt_time14400.txt
	MORE_CONSTANTS=0 time ./verifyAlzette trails/alzette_constb7e15162_rotation8_fullopt_time14400.txt 42 123 456\
		| tee logs/verify_alzette_constb7e15162_rotation8_fullopt_time14400.txt
	MORE_CONSTANTS=0 time ./verifyAlzette trails/alzette_constb7e15162_rotation15_fullopt_time14400.txt 42 123 456\
		| tee logs/verify_alzette_constb7e15162_rotation15_fullopt_time14400.txt
	MORE_CONSTANTS=0 time ./verifyAlzette trails/alzette_constb7e15162_rotation16_fullopt_time14400.txt 42 123 456\
		| tee logs/verify_alzette_constb7e15162_rotation16_fullopt_time14400.txt

	MORE_CONSTANTS=1 time ./verifyAlzette trails/alzette_constb7e15162_rotation1_fullopt_time14400.txt 42 123 456\
		| tee logs/verifyMC_alzette_constb7e15162_rotation1_fullopt_time14400.txt
	MORE_CONSTANTS=1 time ./verifyAlzette trails/alzette_constb7e15162_rotation8_fullopt_time14400.txt 42 123 456\
		| tee logs/verifyMC_alzette_constb7e15162_rotation8_fullopt_time14400.txt
	MORE_CONSTANTS=1 time ./verifyAlzette trails/alzette_constb7e15162_rotation15_fullopt_time14400.txt 42 123 456\
		| tee logs/verifyMC_alzette_constb7e15162_rotation15_fullopt_time14400.txt
	MORE_CONSTANTS=1 time ./verifyAlzette trails/alzette_constb7e15162_rotation16_fullopt_time14400.txt 42 123 456\
		| tee logs/verifyMC_alzette_constb7e15162_rotation16_fullopt_time14400.txt
	
	MORE_CONSTANTS=0 time ./verifyAlzette trails/alzette_constb7e15162_rotation8_fullopt_time14400.txt 45 123 456\
		| tee logs/verify_alzette_constb7e15162_rotation8_fullopt_time14400_data44.txt
	MORE_CONSTANTS=1 time ./verifyAlzette trails/alzette_constb7e15162_rotation8_fullopt_time14400.txt 45 123 456\
		| tee logs/verifyMC_alzette_constb7e15162_rotation8_fullopt_time14400_data44.txt


evolutionFigs:
	python3 process_log_bounds_evolution.py logs/alzette_constb7e15162_rotation6_neqopt_time3600.txt figs/evolution1.png 9999
	python3 process_log_bounds_evolution.py logs/alzette_constb7e15162_rotation6_neqopt_time3600.txt figs/evolution1.pdf 9999
	python3 process_log_bounds_evolution.py logs/alzette_constb7e15162_rotation6_neqopt_time3600.txt figs/evolution2.png 500
	python3 process_log_bounds_evolution.py logs/alzette_constb7e15162_rotation6_neqopt_time3600.txt figs/evolution2.pdf 500


clean:
	rm -f *.o */*.o
	
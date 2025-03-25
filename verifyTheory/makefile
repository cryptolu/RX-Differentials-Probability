CXX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++17 -fopenmp -Wno-deprecated-copy

all: verify1 verify2

Verify-thm-1.c: AL16_patch.diff
	# to avoid license issues, download and patch
	rm -f Verify-thm-1.c Rotational_Cryptanalysis_in_the_Presence_of_Constants.zip
	wget 'http://homes.esat.kuleuven.be/~tashur/Rotational_Cryptanalysis_in_the_Presence_of_Constants.zip'
	unzip Rotational_Cryptanalysis_in_the_Presence_of_Constants.zip Verify-thm-1.c
	patch Verify-thm-1.c AL16_patch.diff

verify1: Verify-thm-1.c
	gcc Verify-thm-1.c -lm -O3 -o verify-thm-1
	./verify-thm-1
	echo "16384 = 2^14 vs 2^15 (2^32 * 2^-17 prob.) claimed by previous theorem"
	echo
	gcc Verify-thm-1.c -lm -O3 -DMISMATCH -o verify-thm-1-mismatch
	./verify-thm-1-mismatch
	echo "1523712 = 2^20.53915881110803 impossible by previous theorem"
	echo

verify2-fast: verify_theory.py
	echo Verifying up to n=8
	# note: can replace pypy3 with python3 but it is 30x slower
	pypy3 verify_theory.py 8

verify2-long: verify_theory.py
	echo Verifying up to n=12
	# note: can replace pypy3 with python3 but it is 30x slower
	pypy3 verify_theory.py 12


verifyTheory: verifyTheory.cpp rxdp.cpp rxdp.hpp common.cpp common.hpp
	$(CXX) $(CXXFLAGS) -o verifyTheory verifyTheory.cpp common.cpp rxdp.cpp

clean:
	rm -f Verify-thm-1.c Rotational_Cryptanalysis_in_the_Presence_of_Constants.zip
	rm -f verify1 verify2-fast verify2-long
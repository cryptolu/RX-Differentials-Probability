CXX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++17 -fopenmp -Wno-deprecated-copy
# CXXFLAGS = -ggdb -O3 -march=native -Wall -Wextra -std=c++17 -fopenmp -Wno-deprecated-copy -g -fsanitize=address -fno-omit-frame-pointer

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

verifyMalzette: verifyMalzette.o
	$(CXX) $(CXXFLAGS) -o verifyMalzette verifyMalzette.o

findInput: findInput.o
	$(CXX) $(CXXFLAGS) -o findInput findInput.o

clean :
	rm -f *.o

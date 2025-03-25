#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <iomanip>

#include "gurobi_c++.h"
#include "common.hpp"

struct Trail {
	std::string cipher;

	uint wordSize;
	uint rotationOffset;
	
	uint nbStateWords;
	uint nbKeyWords;
	
	uint sumNeq;
	double sumLog;

	std::vector<std::vector<uint64_t>> stateWords;
	std::vector<std::vector<uint64_t>> keyWords;
	std::vector<std::vector<uint64_t>> modAddWords;

	void write_to_file(char *filename);
	void read_from_file(char *filename);
};


void setFilename(char * _filename);
void setCipherName(char * cipher);
std::pair<uint, std::vector<std::array<uint64_t,2>>> readKnownTrail(char *_filename);

class CustomCallback: public GRBCallback{

	public:
		uint wordSize;
		uint rotationOffset;
		std::vector<std::array<std::vector<GRBVar>,3>> modAddVars; //list of triplets of list of variables representing the different mod add differentials
		std::vector<std::vector<GRBVar>> stateVars; //List of state variables, listed by round
		GRBVar sumneq;
		std::vector<std::vector<GRBVar>> keyVars; //List of key variables, listed by round, optional

		CustomCallback(uint const wordSize,
					   uint const rotationOffset,
					   std::vector<std::array<std::vector<GRBVar>,3>> const & modAddVars,
					   std::vector<std::vector<GRBVar>> const & stateVars,
					   GRBVar & sumneq,
					   std::vector<std::vector<GRBVar>> const & keyVars = std::vector<std::vector<GRBVar>>());

	protected:
    	void callback();
};
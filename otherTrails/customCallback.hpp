#ifndef H_CUSTOMCALLBACK
#define H_CUSTOMCALLBACK

#include <vector>
#include <array>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <iomanip>

#include "gurobi_c++.h"
#include "common.hpp"

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

#endif
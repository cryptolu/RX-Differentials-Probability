#ifndef H_CUSTOMCALLBACK
#define H_CUSTOMCALLBACK

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <string>
#include <set>

#include "gurobi_c++.h"

class CustomCallback: public GRBCallback{

	public:
		std::vector<GRBVar> keyvar;
		std::vector<unsigned int> cutIndex;
		uint64_t ctrKey;
		std::set<std::vector<uint>> foundKey;

		CustomCallback(std::vector<GRBVar> const & xkeyvar,
					   std::vector<unsigned int> const & xcutIndex);

	protected:
    	void callback();
};

#endif
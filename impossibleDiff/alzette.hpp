#ifndef ALZETTE_H
#define ALZETTE_H

#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <utility>
#include <iostream>
#include <iomanip>
#include <chrono>

#include <omp.h>

#include "gurobi_c++.h"
#include "common.hpp"
#include "MILP_common.hpp"
#include "customCallback.hpp"

//defaultParameters = {31,24,17,17,0,31,24,16};
GRBModel getAlzetteModel(std::vector<uint> const & params,
						 uint64_t const cst,
						 uint const gamma,
						 GRBEnv & env);

void searchWeight1RXDiffAlzette(std::vector<uint> const & params,
								uint64_t const cst);

bool existTrailAlzette(std::vector<uint> const & params,
					   uint64_t const cst,
					   uint const gamma,
					   std::vector<uint> const & inputDiff,
					   std::vector<uint> const & outputDiff,
					   GRBEnv & env);

void checkTruncatedDiffAlzette(std::vector<uint> const & params,
							   uint64_t const cst,
							   uint const gamma,
							   std::vector<uint> const & inputIndex,
							   std::vector<uint> const & outputIndex);

#endif
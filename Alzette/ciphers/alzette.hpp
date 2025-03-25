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
#include <fstream>

#include <omp.h>

#include <assert.h>

#include "gurobi_c++.h"

#include "../common/common.hpp"
#include "../common/rxdp.hpp"
#include "../common/MILP_common.hpp"
#include "../common/customCallback.hpp"

//defaultParameters = {31,24,17,17,0,31,24,16};
GRBModel getAlzetteModel(std::vector<uint> const & params,
						 uint64_t const cst,
						 uint const gamma,
						 GRBEnv & env);
GRBModel getAlzetteModelWithProba(std::vector<uint> const & params,
								  uint64_t const cst,
								  uint const gamma,
								  GRBEnv & env);

std::vector<GRBVar> 
addModAddProbaAlzette(GRBModel & model,
					  std::vector<uint> const & params,
					  uint const gamma,
					  std::vector<std::vector<GRBVar>> & x,
					  std::vector<std::vector<GRBVar>> & y,
					  std::vector<std::vector<GRBVar>> & z);

Trail
searchBestTrailAlzette(std::vector<uint> const & params,
				  uint64_t const cst,
				  uint const gamma,
				  GRBEnv & env,
				  bool const withExactProbability,
				  std::vector<uint64_t> const & inputDiff = std::vector<uint64_t>(),
				  std::vector<uint64_t> const & outputDiff = std::vector<uint64_t>());
//Overload version when no knownTrail is provided
Trail
searchBestTrailAlzette(std::vector<uint> const & params,
				  uint64_t const cst,
				  uint const gamma,
				  GRBEnv & env,
				  bool const withExactProbability,
				  const Trail & knownTrail,
				  std::vector<uint64_t> const & inputDiff = std::vector<uint64_t>(),
				  std::vector<uint64_t> const & outputDiff = std::vector<uint64_t>());
/*
Search for the best RX-diff trail with RX rotation gamma
params should be the rotation values  (default {31,24,17,17,0,31,24,16}), number of rounds is deduced for this as params.size()/2
withExactProbability allows to control if we add probability constraints and optimization. False means we just optimize on the number of neq variables
knownTrail can be provided as a pair (W,trail) where W is a bound on the number of neq variables and trail is a starting solution for the model
inputDiff and outputDiff can be left empty to allow any input/output differential, or given to fix the input/output differential respectively, as {x,y} for the difference on the left,right branch. Won't work for n > 64 since differences are given as a uint64_t for each word (but shouldn't really be a problem).
*/

std::pair<uint32_t,uint32_t> 
alzetteCore(std::vector<uint> const & params,
			uint32_t const cst,
			uint32_t const inx,
			uint32_t const iny);
//Alzette core evaluation

GRBModel
getAlzetteModelPlaintextPair(std::vector<uint> const & params,
							 uint64_t const cst,
							 uint const gamma,
							 GRBEnv & env,
							 const Trail & trail);
//Create a MILP model for Alzette to find a plaintext pair matching a given trail

std::pair<bool, std::array<uint64_t, 4>>
findPlaintextForTrailAlzette(std::vector<uint> const & params,
				  uint64_t const cst,
				  uint const gamma,
				  GRBEnv & env,
				  const Trail & trail);
//Return a pair {b, p}, b = true if a plaintext pair was found (and then p contains the pair as x0,y0,x1,y1)


#endif
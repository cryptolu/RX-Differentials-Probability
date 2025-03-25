#ifndef SPECK_H
#define SPECK_H

#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <utility>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <map>

#include <omp.h>

#include "gurobi_c++.h"
#include "common.hpp"
#include "MILP_common.hpp"
#include "customCallback.hpp"


std::array<std::vector<std::vector<GRBVar>>,2> 
addSpeckKSConstr(GRBModel & model,
				 uint const n,
				 uint const m,
				 uint const nbRound,
				 uint const alpha,
				 uint const beta,
				 std::string const & suffix = "");
//Add constraints for the key schedule of Speck
//These are constraints *in values* not in differences, so constants have an impact
//Return the master key and round key variables

GRBModel
getSpeckModel(uint const n,
			  uint const m, 
			  uint const nbRound,
			  uint const alpha, 
			  uint const beta, 
			  uint const gamma,
			  GRBEnv & env);
//Create a model for single-key RX-diff for Speck 2n/mn over nbRounds rounds, using alpha/beta for the round function. env is only here for better console login from gurobi
//Returns the model

std::vector<GRBVar> 
addModAddProbaSpeck(GRBModel & model,
					uint const alpha,
					uint const gamma,
					std::vector<std::vector<GRBVar>> & x,
					std::vector<std::vector<GRBVar>> & y,
					std::vector<std::vector<GRBVar>> & z);
//Add mod add probability constraints and variables
//Return a vector containing the weight variables

std::tuple<double,uint, std::vector<std::array<uint64_t,3>>>
searchBestTrailSpeck(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				bool const withExactProbability,
				std::map<uint,uint> const & lowerBoundConsecutiveRounds,
				std::vector<uint64_t> const & inputDiff = std::vector<uint64_t>(),
				std::vector<uint64_t> const & outputDiff = std::vector<uint64_t>());
//Overload version when no knownTrail is provided
std::tuple<double,uint, std::vector<std::array<uint64_t,3>>>
searchBestTrailSpeck(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				bool const withExactProbability,
				std::map<uint,uint> const & lowerBoundConsecutiveRounds,
				std::pair<uint, std::vector<std::array<uint64_t,3>>> const & knownTrail,
				std::vector<uint64_t> const & inputDiff = std::vector<uint64_t>(),
				std::vector<uint64_t> const & outputDiff = std::vector<uint64_t>());
/*
Search for the best RX-diff trail with RX rotation gamma
params should be the Speck parameters, in order: n,m,nbRound,alpha,beta
knownTrail can be provided as a pair (W,trail) where W is a bound on the number of neq variables and trail is a starting solution for the model
withExactProbability allows to control if we add probability constraints and optimization
False means we just optimize on the number of neq variables
lowerBoundConsecutiveRounds provides lower bounds on the number of Neq vars for consecutive rounds
e.g. if lowerBoundConsecutiveRounds[2] = 10, then we add constraints so that for any 2 consecutive rounds, there are at least 10 neq vars set to 1
knownTrail can be used to help the solver with a starting point. The first element of the pair is an upper bound on the number of neq vars, and the second element of the pair is the trail given as triplets {x,y,k} for the difference on x,y and the value of k at a given round
inputDiff and outputDiff can be left empty to allow any input/output differential, or given to fix the input/output differential respectively, as {x,y} for the difference on the left,right branch. Won't work for n > 64 since differences are given as a uint64_t for each word (but shouldn't really be a problem).
*/

std::pair<bool, std::array<uint64_t, 4>>
findPlaintextForTrailSpeck(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				std::vector<std::array<uint64_t,3>> const & trail);
/*
//Return a pair {b, p}, b = true if a plaintext pair was found (and then p contains the pair as x0,y0,x1,y1)
*/

std::vector<uint64_t>
speckKS(uint const n,
		uint const nbRound,
		uint const alpha,
		uint const beta,
		std::vector<uint64_t> const & mk);

std::array<uint64_t,2>
speckEncrypt(uint const n,
			 uint const nbRound,
			 uint const alpha,
			 uint const beta,
			 uint64_t x,
			 uint64_t y,
			 std::vector<uint64_t> const & k);

std::array<std::vector<std::vector<GRBVar>>,4> 
addSpeckKSRXDConstr(GRBModel & model,
				 uint const n,
				 uint const m,
				 uint const nbRound,
				 uint const alpha,
				 uint const beta,
				 uint const gamma);
//Add constraints for the key schedule of Speck for RX-differentials
//Return the master key, round key, l and tmpKS variables

std::array<std::vector<std::vector<GRBVar>>,3> 
addSpeckRXDConstrRelatedKey(GRBModel & model,
				  uint const n,
				  uint const nbRound,
				  uint const alpha, 
				  uint const beta, 
				  uint const gamma,
				  std::vector<std::vector<GRBVar>> & k);
//Add constraints for RX-differential propagation in Speck in related key
//Returns a triplet for containing the variables [x,y,z]

GRBModel
getSpeckRelatedKeyModel(uint const n,
			  uint const m, 
			  uint const nbRound,
			  uint const alpha, 
			  uint const beta, 
			  uint const gamma,
			  GRBEnv & env);
//Create a model for related-key RX-diff for Speck 2n/mn over nbRounds rounds, using alpha/beta for the round function. env is only here for better console login from gurobi
//Returns the model

std::pair<std::vector<GRBVar>, std::vector<GRBVar>>
addModAddProbaSpeckRelatedKey(GRBModel & model,
							  uint const alpha,
							  uint const gamma,
							  std::vector<std::vector<GRBVar>> & x,
							  std::vector<std::vector<GRBVar>> & y,
							  std::vector<std::vector<GRBVar>> & z,
							  std::vector<std::vector<GRBVar>> & k,
							  std::vector<std::vector<GRBVar>> & l,
							  std::vector<std::vector<GRBVar>> & tmpKS);
//Add mod add probabilities constraints and variables for the related key model
//Returns a pair of vector of variables, .first the weight variables for the state, .second for the weight variables of the key schedule

std::tuple<double,uint, std::vector<std::array<uint64_t,3>>> 
searchBestTrailSpeckRelatedKey(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				bool const withExactProbability,
				std::map<uint,uint> const & lowerBoundConsecutiveRounds,
				std::vector<uint64_t> const & inputDiff = std::vector<uint64_t>(),
				std::vector<uint64_t> const & outputDiff = std::vector<uint64_t>());

std::tuple<double,uint, std::vector<std::array<uint64_t,3>>> 
searchBestTrailSpeckRelatedKey(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				bool const withExactProbability,
				std::map<uint,uint> const & lowerBoundConsecutiveRounds,
				std::pair<uint, std::vector<std::array<uint64_t,3>>> const & knownTrail,
				std::vector<uint64_t> const & inputDiff = std::vector<uint64_t>(),
				std::vector<uint64_t> const & outputDiff = std::vector<uint64_t>());
//Same as searchBestTrailSpeck but for related key
//the returned trail/the one used as input will contain the key RXdiff trail instead of key value

GRBModel
getSpeckModelPlaintextPairRelatedKey(uint const n,
						   uint const m, 
						   uint const nbRound,
						   uint const alpha, 
						   uint const beta, 
						   uint const gamma,
						   GRBEnv & env,
						   std::vector<std::array<uint64_t,3>> const & trail);
//Create a MILP model for Speck to find a plaintext pair matching a given trail in related key

bool
findPlaintextForTrailSpeckRelatedKey(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				std::vector<std::array<uint64_t,3>> const & trail);

#endif
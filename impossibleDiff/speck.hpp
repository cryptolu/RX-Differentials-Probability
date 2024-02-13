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

#include <omp.h>

#include "gurobi_c++.h"
#include "common.hpp"
#include "MILP_common.hpp"
#include "customCallback.hpp"

void addSpeckKSConstr(GRBModel & model,
					  std::vector<std::vector<GRBVar>> & mk,
					  std::vector<std::vector<GRBVar>> & k,
					  uint const alpha,
					  uint const beta);
//Add constraints for the key schedule of Speck
//mk are the variables for the master key (split in m vectors)
//n (word size) is obtained from mk[0].size()
//These are constraints *in values* not in differences, so constants have an impact

GRBModel
getSpeckModel(uint const n,
			  uint const m, 
			  uint const nbRound,
			  uint const alpha, 
			  uint const beta, 
			  uint const gamma,
			  GRBEnv & env);
//Create a model for single-key RX-diff for Speck 2n/mn over nbRounds rounds, using alpha/beta for the round function. env is only here for better console login from gurobi
GRBModel
getSpeckModelNoKey(uint const n,
				   uint const nbRound,
				   uint const alpha, 
				   uint const beta, 
				   uint const gamma,
				   GRBEnv & env);

void searchWeight1RXDiffSpeck(uint const n,
							  uint const m,
							  uint const nbRound,
							  std::vector<uint> const & cutIndex);

void searchWeight1RXDiffSpeckNoKey(uint const n,
								   uint const nbRound);

bool existTrailSpeckNoKey(uint const n,
						  uint const nbRound,
						  uint const gamma,
						  std::vector<uint> const & inputDiff,
						  std::vector<uint> const & outputDiff,
						  GRBEnv & env);



void addSpeckKSRelatedKeyConstr(GRBModel & model,
								std::vector<std::vector<GRBVar>> & mk,
								std::vector<std::vector<GRBVar>> & k,
								uint const alpha,
								uint const beta,
								uint const gamma);

GRBModel
getSpeckRelatedKeyModel(uint const n,
						uint const m, 
						uint const nbRound,
						uint const alpha, 
						uint const beta, 
						uint const gamma,
						GRBEnv & env);

void searchWeight1RXDiffSpeckRelatedKey(uint const n,
										uint const m,
										uint const nbRound);


bool existTrailSpeckRelatedKey(uint const n,
							   uint const m,
							   uint const nbRound,
							   uint const gamma,
							   std::vector<uint> const & inputDiff,
							   std::vector<uint> const & outputDiff,
							   std::vector<uint> const & keyDiff,
							   GRBEnv & env);

void checkTruncatedDiffSpeckRelatedKey(uint const n,
									   uint const m,
									   uint const nbRound,
									   uint const gamma,
									   std::vector<uint> const & inputIndex,
									   std::vector<uint> const & outputIndex,
									   std::vector<uint> const & keyIndex);
#endif
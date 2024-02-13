#ifndef XTEA_H
#define XTEA_H

#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <utility>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "gurobi_c++.h"
#include "common.hpp"
#include "MILP_common.hpp"
#include "customCallback.hpp"

void addXteaKSConstr(GRBModel & model,
					std::vector<std::vector<GRBVar>> & mk,
					std::vector<std::vector<GRBVar>> & k);
//Add constraints for the XTEA key schedule
//mk are the variables for the master key (split in m vectors of 32 vars)
//k are the round key variables, split in r*32 vars
//These are constraints *in values* not in differences, so constants have an impact

GRBModel getXteaModel(uint const nbRound,
					  uint const gamma,
					  GRBEnv & env);

GRBModel getXteaModelNoKEy(uint const nbRound,
						   uint const gamma,
						   GRBEnv & env);

void searchWeight1RXDiffXtea(uint const nbRound,
							 std::vector<uint> const & cutIndex);
void searchWeight1RXDiffXteaNoKey(uint const nbRound);
#endif
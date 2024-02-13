#ifndef LEA_H
#define LEA_H

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

void addLeaKSConstr(GRBModel & model,
					std::vector<std::vector<GRBVar>> & mk,
					std::vector<std::vector<std::vector<GRBVar>>> & k);
//Add constraints for the key schedule of LEA
//mk are the variables for the master key (split in m vectors of 32 vars)
//k are the round key variables, split in r*6*32 vars
//These are constraints *in values* not in differences, so constants have an impact

GRBModel getLeaModel(uint const keySize,
					 uint const nbRound,
					 uint const gamma,
					 GRBEnv & env);
GRBModel getLeaModelNoKey(uint const nbRound,
						  uint const gamma,
						  GRBEnv & env);

void searchWeight1RXDiffLea(uint const keySize,
							uint const nbRound,
							std::vector<uint> const & cutIndex);
void searchWeight1RXDiffLeaNoKey(uint const nbRound);

#endif
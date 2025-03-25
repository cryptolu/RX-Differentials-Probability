#ifndef SALSA_H
#define SALSA_H

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

std::vector<GRBVar>
addSalsaOpConstr(GRBModel & model,
				 std::vector<GRBVar> & a,
				 std::vector<GRBVar> & b,
				 std::vector<GRBVar> & c,
				 std::vector<GRBVar> & out,
				 uint const rot,
				 uint const gamma,
				 std::string const & prefix);
/*
Add the constraints to modelize RX-differentials of out = ((a + b) <<< rot) ^ c
Create temp variables named prefix+"_"+i for i in [0,31]
Return the vector of temp vars
*/

std::vector<std::vector<GRBVar>>
addSalsaQuarterRoundConstr(GRBModel & model,
						   uint const gamma,
						   std::array<std::vector<GRBVar>,4> & in,
						   std::array<std::vector<GRBVar>,4> & out,
						   std::string const & suffix);
/*
Add constraints between in and out to modelize the RXDiff propagation from in to out
in and out should be of size (4,32)
suffix is used for unique variable names
Return an array of 4 vectors of variables used to modelize the mod add
Quarter round is defined as 
input {a,b,c,d}
b ^= (a + d) <<< 7;
c ^= (b + a) <<< 9;
d ^= (c + b) <<< 13;
a ^= (d + c) <<< 18;
output {a,b,c,d}
*/

std::pair<std::vector<std::vector<std::vector<GRBVar>>>, 
		  std::vector<std::vector<std::vector<std::vector<GRBVar>>>>>
addSalsaRXDConstr(GRBModel & model,
				  uint const nbRound,
				  uint const startingRound,
				  uint const gamma);
/*
Add constraints for RX-differential propagation for Salsa
Return the variables {x,z} used in the modelization
*/

void addNeqConstrSingle(GRBModel & model,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						std::vector<GRBVar> & z,
						std::vector<GRBVar> & neq,
						uint const gamma);
//Add constraints neq[i] = NotAllEqual(x[i],y[i],z[i]) for all i != (gamma-1)

void 
addSalsaNeqConstr(GRBModel & model,
				  uint const nbRound,
				  uint const startingRound,
				  uint const gamma,
				  std::vector<std::vector<std::vector<GRBVar>>> & x,
				  std::vector<std::vector<std::vector<std::vector<GRBVar>>>> & z,
				  std::vector<std::vector<std::vector<std::vector<GRBVar>>>> & neq);
//Link neq vars to state vars

GRBModel
getSalsaModel(uint const nbRound,
			  uint const startingRound,
			  uint const gamma,
			  GRBEnv & env);
/*
Generate a model of nbRound rounds of Salsa, starting from startingRound
startingRound can be used to adjust whether the first round should be even or odd
*/

std::vector<std::vector<std::vector<GRBVar>>>
addModAddProbaSalsa(GRBModel & model,
					uint const nbRound,
					uint const startingRound,
					uint const gamma,
					std::vector<std::vector<std::vector<GRBVar>>> & x,
					std::vector<std::vector<std::vector<std::vector<GRBVar>>>> & z);

std::tuple<double,uint, std::vector<std::array<uint64_t,16>>> 
searchBestTrailSalsa(uint const nbRound,
					 uint const startingRound,
					 uint const gamma,
					 GRBEnv & env,
					 bool const withExactProbability,
					 std::map<uint,uint> const & lowerBoundConsecutiveRounds,
					 std::vector<uint64_t> const & inputDiff = std::vector<uint64_t>(),
					 std::vector<uint64_t> const & outputDiff = std::vector<uint64_t>());

std::tuple<double,uint, std::vector<std::array<uint64_t,16>>> 
searchBestTrailSalsa(uint const nbRound,
					 uint const startingRound,
					 uint const gamma,
					 GRBEnv & env,
					 bool const withExactProbability,
					 std::map<uint,uint> const & lowerBoundConsecutiveRounds,
					 std::pair<uint, std::vector<std::array<uint64_t,16>>> const & knownTrail,
					 std::vector<uint64_t> const & inputDiff = std::vector<uint64_t>(),
					 std::vector<uint64_t> const & outputDiff = std::vector<uint64_t>());
/*
Search for the best RX-diff trail with RX rotation gamma
withExactProbability allows to control if we add probability constraints and optimization. False means we just optimize on the number of neq variables
knownTrail can be provided as a pair (W,trail) where W is a bound on the number of neq variables and trail is a starting solution for the model
inputDiff and outputDiff can be left empty to allow any input/output differential, or given to fix the input/output differential respectively
*/

bool
findPlaintextForTrailSalsa(uint const nbRound, 
						   uint const startingRound,
						   uint const gamma,
						   GRBEnv & env,
						   std::vector<std::array<uint64_t,16>> const & trail,
						   bool const printOutput = true);


class SalsaCallback: public GRBCallback{

	public:
		uint wordSize;
		uint rotationOffset;
		std::vector<std::array<std::vector<GRBVar>,3>> modAddVars; //list of triplets of list of variables representing the different mod add differentials
		std::vector<std::vector<GRBVar>> stateVars; //List of state variables, listed by round
		GRBVar sumneq;
		uint nbRound;
		uint startingRound;
		GRBEnv env;
		uint64_t nbSolFound;

		SalsaCallback(uint const wordSize,
					  uint const rotationOffset,
					  std::vector<std::array<std::vector<GRBVar>,3>> const & modAddVars,
					  std::vector<std::vector<GRBVar>> const & stateVars,
					  GRBVar & sumneq,
					  uint const nbRound,
					  uint const startingRound,
					  GRBEnv & env);

	protected:
    	void callback();
};


std::tuple<double,uint, std::vector<std::array<uint64_t,4>>> 
searchBestTrailSalsaQR(int const indexCst,
					   uint const gamma,
					   GRBEnv & env,
					   bool const withExactProbability);

std::tuple<double,uint, std::vector<std::array<uint64_t,4>>> 
searchBestTrailSalsaQR(int const indexCst,
					   uint const gamma,
					   GRBEnv & env,
					   bool const withExactProbability,
					   std::pair<uint, std::vector<std::array<uint64_t,4>>> const & knownTrail);
//Search for the best RX-diff trail over one quarter round 
//indexCst = -1 means we don't use the RXD of constants as constraints
//otherwise, fix the word of index indexCst to its corresponding RXD


bool
findPlaintextForTrailSalsaQR(int const indexCst,
							 uint const gamma,
							 GRBEnv & env,
							 std::vector<std::array<uint64_t,4>> const & trail,
							 bool const printOutput = true);



class SalsaQRCallback: public GRBCallback{

	public:
		uint wordSize;
		uint rotationOffset;
		std::vector<std::array<std::vector<GRBVar>,3>> modAddVars; //list of triplets of list of variables representing the different mod add differentials
		std::vector<std::vector<GRBVar>> stateVars; //List of state variables, listed by round
		GRBVar sumneq;
		uint nbRound;
		GRBEnv env;
		uint64_t nbSolFound;
		int indexCst;

		SalsaQRCallback(uint const wordSize,
					  uint const rotationOffset,
					  std::vector<std::array<std::vector<GRBVar>,3>> const & modAddVars,
					  std::vector<std::vector<GRBVar>> const & stateVars,
					  GRBVar & sumneq,
					  uint const nbRound,
					  GRBEnv & env,
					  int const indexCst);

	protected:
    	void callback();
};


#endif
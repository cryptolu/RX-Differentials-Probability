#ifndef MILP_COMMON_H
#define MILP_COMMON_H

#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <map>

#include "gurobi_c++.h"
#include "common.hpp"

using uint = unsigned int;

/*---- XOR Constraints ----
Separate implementations for 2 and 3 variables
Dummy versions for 2 & 3 are there for completness but should probably be avoided
Anything with >= 4 variables will use the dummy version
*/
//-- 2 input variables --
void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & z,
				  uint const c=0);
//Add constraint z = x XOR y XOR c with x,y,z binaries, c a constant (0 or 1)
void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & z,
				  std::string const & dumName);
//Add constraint z = x XOR y with x,y,z binaries
//Dummy variable version, probably better not to use
void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & z,
				  uint const c,
				  std::string const & dumName);
//Add constraint z = x XOR y XOR c with x,y,z binaries, c constant
//Dummy variable version, probably better not to use

//-- 3 input variables --
void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & w,
				  GRBVar & z,
				  uint const c=0);
//Add constraint z = x XOR y XOR w XOR c with x,y,w,z binaries, c constant
void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & w,
				  GRBVar & z,
				  std::string const & dumName);
//Add constraint z = x XOR y XOR w with x,y,w,z binaries
//Dummy variable version, probably better not to use
void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & w,
				  GRBVar & z,
				  uint const c,
				  std::string const & dumName);
//Add constraint z = x XOR y XOR w XOR c with x,y,w,z binaries, c constant
//Dummy variable version, probably better not to use

//-- n input variables --
void addXORConstr(GRBModel & m,
				  std::vector<GRBVar> & vars,
				  GRBVar & z,
				  std::string const & dumName);
//Add constraint z = vars[0] XOR vars[1] XOR ... XOR vars[n-1]
//Arbitrary length version, always use dummy variable if >= 4 input variables
//if 2 or 3 input variables, uses the binary modelization for now as it's likely better
void addXORConstr(GRBModel & m,
				  std::vector<GRBVar> & vars,
				  GRBVar & z,
				  uint const c,
				  std::string const & dumName);
//Add constraint z = vars[0] XOR vars[1] XOR ... XOR vars[n-1] XOR c
//Arbitrary length version, always use dummy variable if >= 4 input variables
//if 2 or 3 input variables, uses the binary modelization for now as it's likely better

//--RX-diff XOR with a constant --
//RX-differences behave differently when XOR'd with a constant, which could allow us to consider trails in the single key model by having the key schedule modelized in values, and the state as differences
void addRXCstXORConstr(GRBModel & m,
					   std::vector<GRBVar> & vars,
					   std::vector<GRBVar> & cst,
					   std::vector<GRBVar> & res,
					   uint const s);
//Add a constraint for the RX-diff propagation through res = vars XOR cst where cst is seen as a constant but still a gurobi var (e.g. round key with key-schedule modelized) with RX shift s

void addRXCstXORConstr(GRBModel & m,
					   std::vector<GRBVar> & vars,
					   std::vector<uint> const & cst,
					   std::vector<GRBVar> & res,
					   uint const s);
//Add a constraint for the RX-diff propagation through res = vars XOR cst where cst is a constant (given as a fixed binary vector) with RX shift s


/*---- Simon-like round function ----
This is used to modelize differences through the function f(x) = (x << a) & (x << b) ^ (x << c)
It's agnostic to the shift in an RX differential, so no need to consider it there
*/
void addSimonCoreConstr(GRBModel & m,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						int a,
						int b,
						int const c,
						std::string const & prefix);
//Add constraints to modelize the differential transition x -> y through the function
//f(x) = (x << a) & (x << b) ^ (x << c)
//x and y are assumed to be of the same size
//Requires the introduction of variables with varname starting by prefix

void addSimonCoreValConstr(GRBModel & model,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						int a,
						int b,
						int const c);
//Add constraints to modelize the transition **in value** y = f(x)
//This version does a direct modelization without intermediate variables


/*---- Mod Add Constraints ----
Constraints related to modular addition
*/
void addModAddValueConstr(GRBModel & model, 
						  std::vector<GRBVar> & x,
						  std::vector<GRBVar> & y,
						  std::vector<GRBVar> & z,
						  std::string const & prefix);
//Constraints to represent z = x + y **in value**
//prefix is used to create carry variables

void addModAddConstValueConstr(GRBModel & model, 
							   std::vector<GRBVar> & x,
							   uint64_t const cst,
							   std::vector<GRBVar> & z,
							   std::string const & prefix);
//Constraints to represent z = x + y **in value** where y is a known constant
//prefix is used to create carry variables
//Limited to constants over 64 bits for now

void fConstraintRXDiffModAdd(GRBModel & model,
							 GRBVar & x1,
							 GRBVar & x2,
							 GRBVar & x3,
							 GRBVar & x4,
							 GRBVar & x5,
							 GRBVar & x6);
//f-constraint used in the modelization of the RX-diff propagation for mod add

void addModAddRXDiffConstr(GRBModel & model,
						   std::vector<GRBVar> & a,
						   std::vector<GRBVar> & b,
						   std::vector<GRBVar> & d,
						   uint const gamma);
//Add constraints to modelize the RX-diff propagation for mod add
// (a,gamma),(b,gamma) -> (d,gamma)


/*
---- Normal (non-cyclic) shift RX-Diff constraints ----
*/

void addSHLRXDiffConstr(GRBModel & model,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						int const alpha,
						int const gamma);
//Add constraints to modelize the RX-Diff propagation for y = x << alpha


void addSHRRXDiffConstr(GRBModel & model,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						int const beta,
						int const gamma);
//Add constraints to modelize the RX-Diff propagation for y = x << beta


//--------------------------------------------
//---- Constraints for RXDP probabilities ----
//--------------------------------------------

void addAllEqualBoolConstr(GRBModel & model,
						   std::vector<GRBVar> & x,
						   GRBVar & eq);
//Add a constraints so that eq == 1 iif all variables in x have the same value
void addNotAllEqualBoolConstr(GRBModel & model,
							  std::vector<GRBVar> & x,
							  GRBVar & eq);
//Add a constraints so that eq == 0 iif all variables in x have the same value

void modelMux(GRBModel & model,
			  GRBVar & dst,
			  GRBVar & flag,
			  GRBLinExpr & expr0,
			  GRBLinExpr & expr1,
			  double const M);
void modelMux(GRBModel & model,
			  GRBVar & dst,
			  GRBVar & flag,
			  GRBLinExpr & expr0,
			  GRBLinExpr & expr1);

void addTnConstr(GRBModel & model,
				 std::vector<GRBVar> & alpha,
				 std::vector<GRBVar> & beta,
				 std::vector<GRBVar> & delta,
				 GRBVar & lsbOther,
				 GRBVar & wt,
				 std::string const & suffix);
	/*
	Add constraints to the model to modelize the relation
	wt = -log2(Tn(alpha,beta,delta,lsbOther))
	suffix is used to create intermediary variables with unique names:
	n is obtained from the size of the alpha vector
	*/

void addModAddRXProbaConstr(GRBModel & model,
							std::vector<GRBVar> & alpha,
							std::vector<GRBVar> & beta,
							std::vector<GRBVar> & delta,
							GRBVar & wt,
							int const k,
							std::string const & suffix);
	/*
	Add constraints to modelize the probability of the RX-differential (alpha,beta) -> delta
	through the modular addition, with RX-rotation k
	The result is given by wt = -log2(probability)
	suffix is used to create intermediary variables with unique names
	n is obtained from the size of the alpha vector
	*/

void addVectorEqualValConstr(GRBModel & model,
							 std::vector<GRBVar> & vars,
							 std::vector<uint> const & val);
/*
Add constraints so that vars[i] == val[i]
*/

uint64_t getUintSolutionFromBinVector(std::vector<GRBVar> & vars);
//Given a vector of binary variables as input, return the corresponding integer

std::vector<GRBVar> getArrayVar(GRBModel & model,
								uint const dim,
								std::string const & prefix);
//Get all vars in an array s.t. t[i] is the variable names prefix+i

std::vector<std::vector<GRBVar>> 
get2DArrayVar(GRBModel & model,
			  uint const dim1,
			  uint const dim2,
			  std::string const & prefix);
//Get all vars in a 2D array s.t. t[i][j] is the variable names prefix+i+_+j

std::vector<GRBVar> genArrayBinVar(GRBModel & model,
								   uint const dim,
								   std::string const & prefix);
//Return a vector of binary variables named prefix+i

std::vector<std::vector<GRBVar>> genArray2DBinVar(GRBModel & model,
												  uint const dim1,
												  uint const dim2,
												  std::string const & prefix);
//Return a 2D vector of binary variables named prefix+i+_+j

void addRXDifferentialValueBinding(GRBModel & model,
								   std::vector<GRBVar> & x0,
								   std::vector<GRBVar> & x1,
								   std::vector<GRBVar> & dx,
								   uint const gamma);
//Add the constraints to enforce dx = rotl(x0,gamma) xor x1

void addValueToBinVectorConstr(GRBModel & model,
							   std::vector<GRBVar> & vars,
							   uint64_t val);
//add constraints vars[i] == val[i] (val[i] = i-th bit of val)

void
addMatsuiLikeConstr(GRBModel & model,
					std::map<uint,uint> const & lowerBoundConsecutiveRounds,
					std::vector<GRBVar> & vars,
					GRBVar & sumneq);
//Add matsui-like constraints to derive lower bounds on consecutive rounds
//e.g. if lowerBoundConsecutiveRounds[3] = 5, then any sum of three consecutive variables in vars must be >= 5
//Also add constraints from ZSCH2018 but those might be redundant

#endif
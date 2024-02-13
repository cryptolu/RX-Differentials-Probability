#ifndef MILP_COMMON_H
#define MILP_COMMON_H

#include <string>
#include <vector>
#include <cmath>
#include <numeric>

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

#endif
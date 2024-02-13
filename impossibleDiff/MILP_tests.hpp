#ifndef MILP_TESTS_H
#define MILP_TESTS_H

#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <map>
#include <utility>
#include <random>

#include "gurobi_c++.h"
#include "MILP_common.hpp"
#include "common.hpp"
#include "speck.hpp"
#include "lea.hpp"
#include "xtea.hpp"

using uint = unsigned int;

bool testXOR2(bool const dummyVar, uint const c=0);
bool testXOR3(bool const dummyVar, uint const c=0);
bool testXORn(uint const nvar, uint const c=0);
bool testRXCstXOR(uint const nvar);
bool testRXCstXOR_nonCstVar(uint const nvar);

bool testModAddValue(uint const n);
bool testModAddConstValue(uint const n);

bool testSpeckKS(uint const n);

bool testModAddRXDiff(uint const n,
					  uint const gamma);

bool testSpeckValue(uint const n,
					uint const m);

bool testLeaKS(uint const keySize);

bool testLeaValue(uint const keySize);

std::pair<uint32_t,uint32_t> alzetteCore(std::vector<uint> const & params,
										 uint32_t const cst,
										 uint32_t const inx,
										 uint32_t const iny);
bool testAlzetteValue(std::vector<uint> const & params,
					  uint const cst);

bool testXteaKS(std::vector<uint32_t> const & key);
void xteaEncrypt(unsigned int num_rounds, uint32_t v[2], uint32_t const key[4]);
bool testXteaValue();

bool testSHLRXDiff(int const n,
				   int const alpha,
				   int const gamma);

bool testSHRRXDiff(int const n,
				   int const beta,
				   int const gamma);

bool checkCriteria(uint64_t const a,
				   uint64_t const b,
				   uint64_t const d,
				   uint const n,
				   uint const gamma);

bool testTheoremModAdd(uint const n,
					   uint const gamma);

#endif
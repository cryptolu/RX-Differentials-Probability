#ifndef TESTSHIFTRXDIFF_H
#define TESTSHIFTRXDIFF_H

#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <map>
#include <utility>
#include <random>
#include <sstream>
#include <string>

#include "common.hpp"

//Value tests are done bit by bit to exactly reproduce the formulas, so it is rather slow
//Could be done with masks instead, which is done for RX-diff tests, but those are the ones that really matters anyway
uint getBit(uint64_t const x,
			int const i,
			uint const n);

void checkValueSHL(int const n,
				   int const gamma,
				   int const alpha,
				   uint64_t const x0,
				   uint64_t const x1,
				   uint64_t const deltaOut);
void testValueSHL(uint const n);

void checkValueSHR(int const n,
				   int const gamma,
				   int const beta,
				   uint64_t const x0,
				   uint64_t const x1,
				   uint64_t const deltaOut);
void testValueSHR(uint const n);

bool checkRXDiffSHL(int const n,
					int const gamma,
					int const alpha,
					uint64_t const deltaIn,
					uint64_t const deltaOut);
void testRXDiffSHL(uint const n);

bool checkRXDiffSHR(int const n,
					int const gamma,
					int const beta,
					uint64_t const deltaIn,
					uint64_t const deltaOut);
void testRXDiffSHR(uint const n);

#endif
#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <iostream>
#include <set>
#include <vector>
#include <fstream>
#include <string>

using uint = unsigned int;

uint64_t CSHL(uint64_t const x,
		  uint const s,
		  uint const n);


uint64_t CSHR(uint64_t const x,
		  uint const s,
		  uint const n);

uint64_t shl(uint64_t const x,
			 uint const s,
			 uint const n);

uint64_t shr(uint64_t const x,
			 uint const s,
			 uint const n);

void binprint(uint64_t const x,
			  uint const n,
			  bool const rightLSB = true);

int mod(int const x,
		int const n);

uint64_t binVectorToUint(std::vector<uint> const & v);
uint64_t nextUintSameHammingWeight(uint const val);

std::set<uint64_t> 
genAllVectorHammingDistance(std::vector<std::vector<uint>> const & knownVectors,
							uint const distance);

void writeVectorUint64ToFile(std::vector<uint64_t> const & v,
							 std::string const & filename,
							 std::ios_base::openmode mode = std::ios_base::out);

std::vector<uint64_t> readVectorUint64FromFile(std::string const & filename,
											   std::ios_base::openmode mode = std::ios::in,
											   uint64_t const offset=0);



bool checkCriteria(uint64_t const a,
				   uint64_t const b,
				   uint64_t const d,
				   uint const n,
				   uint const k);

uint64_t Tcounter(uint64_t const alpha,
				  uint64_t const beta,
				  uint64_t const delta,
				  uint const n,
				  uint const v);

uint64_t getRXDiffCount(uint64_t const alpha,
						uint64_t const beta,
						uint64_t const delta,
						uint const n,
						uint const k);

void checkCountFormula(uint const n, uint const k);

#endif
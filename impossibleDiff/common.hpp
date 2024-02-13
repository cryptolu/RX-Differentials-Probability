#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <iostream>

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
			  bool const rightLSB = false);

int mod(int const x,
		int const n);

#endif
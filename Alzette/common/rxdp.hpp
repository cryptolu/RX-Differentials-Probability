#pragma once

#include "common.hpp"

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

uint64_t RXDPconstTcounter_i(uint64_t a, uint64_t b, uint64_t delta, uint64_t Delta, uint64_t C, int n, int i);

uint64_t getRXDiffConstCount(uint64_t a, uint64_t delta, uint64_t Delta, int n, int k);
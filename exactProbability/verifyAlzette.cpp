#include <bits/stdc++.h>
#include <unistd.h>

#include "common/customCallback.hpp"

static inline uint32_t rol(uint32_t x, unsigned int n) {
	return (x<<n) | (x>>(32-n));
}
static inline uint32_t ror(uint32_t x, unsigned int n) {
	return (x>>n) | (x<<(32-n));
}

using namespace std;

bool INDEPENDENT_ROUNDS=false;
bool MORE_CONSTANTS=false;
Trail trail;

int ROUNDS = 4;
int K = 3;
uint32_t CONST;

const uint32_t consts_Alzette[] = {
	0xB7E15162,0xBF715880,0x38B4DA56,0x324E7738,
	0xBB1185EB,0x4F7C7B57,0xCFBFA1C8,0xC2B3293D,
};

const pair<int, int> rotations[] = {
    {31, 24},
    {17, 17},
    {0, 31},
    {24, 16},
};

static inline uint32_t hash32(uint64_t x) {
	x += 0xa40a092dc6d14103ull;
	x *= 0x6353cd1ebf1eea4dull;
	x ^= x >> 17;
	x ^= x >> 37;
	x ^= x >> 23;

	x += 0xb3fd57a50cd6548full;
	x *= 0xf0646bbe961da6afull;
	x ^= x >> 11;
	x ^= x >> 39;
	x ^= x >> 19;
	return x;
}

void Alzette(uint32_t &l, uint32_t &r, uint32_t cst, bool pairParity, uint32_t itr) {
	uint32_t ref;
	if (INDEPENDENT_ROUNDS)
		ref = hash32(itr);
	
	uint32_t negs;
	if (MORE_CONSTANTS)
		negs = hash32(hash32(itr));

	for(int i = 0; i < ROUNDS; i++) {
		l += ror(r, rotations[i].first);
		r ^= ror(l, rotations[i].second);
		l ^= cst;

		if (INDEPENDENT_ROUNDS) {
			// randomize to test round dependency effect
			ref ^= 0xabcd1234;
			ref ^= (~ref) >> 3;
			ref ^= (~ref) << 11;
			l ^= pairParity ? rol(ref, K) : ref;
			ref ^= (~ref) >> 3;
			ref ^= (~ref) << 11;
			r ^= pairParity ? rol(ref, K) : ref;
		}
		if (MORE_CONSTANTS) {
			if (negs & 1) l = ~l;
			if (negs & 2) r = ~r;
			negs >>= 2;
		}
	}
}

bool followsTrail(uint32_t l, uint32_t r, uint32_t ll, uint32_t rr, uint32_t cst) {
	for(int i = 0; i < ROUNDS; i++) {
		if (1) { // test
			assert((rol(l, K) ^ ll) == trail.modAddWords[i][0]);
			uint32_t r1 = ror(r, rotations[i].first);
			uint32_t r2 = ror(rr, rotations[i].first);
			assert((rol(r1, K) ^ r2) == trail.modAddWords[i][1]);
		}
		l += ror(r, rotations[i].first);
		ll += ror(rr, rotations[i].first);

		uint32_t delta = rol(l, K) ^ ll;
		if (delta != trail.modAddWords[i][2]) {
			printf("  failed follow at round %d/%d: %08x != %08lx\n", i+1, ROUNDS, delta, trail.modAddWords[i][2]);
			return false;
		}

		r ^= ror(l, rotations[i].second);
		rr ^= ror(ll, rotations[i].second);
		
		l ^= cst;
		ll ^= cst;
	}
	return true;
}

int main(int argc, char *argv[]) {
	srandom(time(NULL));
	srandom(random() ^ getpid());
	uint64_t seed1 = (random() << 32) ^ random();
	uint64_t seed2 = (random() << 32) ^ random();

	if (argc-1 < 2) {
		printf("Usage: %s <trail_file> <log2(#pairs)> <seed1:hex> <seed2:hex>\n", argv[0]);
		return 0;
	}

	int data_log = 32;

	trail.read_from_file(argv[1]);
	data_log = atoi(argv[2]);

	if (argc-1 >= 4) {
		seed1 = strtoull(argv[3], NULL, 16);
		seed2 = strtoull(argv[4], NULL, 16);
	}

	char *indep_rounds_flag=getenv("INDEPENDENT_ROUNDS");
	if (indep_rounds_flag) {
		INDEPENDENT_ROUNDS = atoi(indep_rounds_flag);
	}
	char *more_constants_flag=getenv("MORE_CONSTANTS");
	if (more_constants_flag) {
		MORE_CONSTANTS = atoi(more_constants_flag);
	}

	printf("===========================\n");
	if (INDEPENDENT_ROUNDS) {
		printf("INDEPENDENT_ROUNDS TESTING!\n");
	}
	else if (MORE_CONSTANTS) {
		printf("MORE CONSTANTS TESTING\n");
	}
	else {
		printf("Normal rounds testing\n");
	}
	printf("===========================\n");

	assert(3 == sscanf(trail.cipher.c_str(), "Alzette_%08x_k%d_nr%d", &CONST, &K, &ROUNDS));

	assert(ROUNDS == 4);
	assert(1 <= K && K <= 16);
	// K = 32-K;

	printf("Testing %s, using 2^%d data, seeds %016lx %016lx\n", trail.cipher.c_str(), data_log, seed1, seed2);
	printf("Trail (state words):\n");
	for (auto &state: trail.stateWords) {
		printf("  ");
		for (auto word: state) {
			printf("%08x ", word);
		}
		printf("\n");
	}

	uint64_t hits = 0;
	uint64_t hits_trail = 0;
	uint64_t n_iters = 1ull << data_log;
	
	ROUNDS = 4;
	auto inputDiff = trail.stateWords[0];
	auto outDiff = trail.stateWords[ROUNDS];

	// randomize human entered seed
	seed1 ^= hash32(seed2);
	seed1 ^= hash32(seed1);
	seed2 ^= hash32(seed1);
	seed2 ^= hash32(seed2);

	#pragma omp parallel for
	for (uint64_t itr = 0; itr < n_iters; itr++) {
		uint32_t l0, r0, l, r, ll, rr, ll0, rr0;
		// note: random() slows openmp a lot!
		// l0 = random() ^ (random() << 16);
		// r0 = random() ^ (random() << 16);

		// use custom random if openmp
		l = l0 = hash32(seed1 ^ itr);
		r = r0 = hash32(seed2 ^ itr);
		Alzette(l, r, CONST, 0, itr);

		ll = ll0 = rol(l0, K) ^ inputDiff[0];
		rr = rr0 = rol(r0, K) ^ inputDiff[1];
		Alzette(ll, rr, CONST, 1, itr);

		if (ll != (rol(l, K) ^ outDiff[0]))
			continue;
		if (rr != (rol(r, K) ^ outDiff[1]))
			continue;
		
		#pragma omp critical
		{
			hits++;
			// bool follows = followsTrail(l0, r0, ll0, rr0, CONST);
			bool follows = 1;
			hits_trail += follows;
			printf("Found match #%lu %08x:%08x (follows trail? %d)\n",
				hits, l0, r0, follows
			);
			fflush(stdout);
		}
	}
	printf("%s: %lu/%lu = 2^%.2f experimental differential probability\n",
		trail.cipher.c_str(),
		hits, n_iters, (log(hits)-log(n_iters))/log(2));
	printf("%s: %lu/%lu/%lu = 2^%.2f experimental trail probability\n",
		trail.cipher.c_str(),
		hits_trail, hits, n_iters, (log(hits_trail)-log(n_iters))/log(2));
	printf("\n");
	fflush(stdout);
	return 0;
}
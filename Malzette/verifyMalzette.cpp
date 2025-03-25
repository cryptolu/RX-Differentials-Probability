#include <bits/stdc++.h>
#include <unistd.h>

static inline uint32_t rol(uint32_t x, unsigned int n) {
	return (x<<n) | (x>>(32-n));
}
static inline uint32_t ror(uint32_t x, unsigned int n) {
	return (x>>n) | (x<<(32-n));
}

using namespace std;

const int ROUNDS = 12;
const int K = 3;

const pair<uint32_t, uint32_t> consts_Malzette1[] = {
    {0x00000000, 0x4e381c1c},
    {0x2aaaaaaa, 0x36dbe492},
    {0x7fffffff, 0x1236db6c},
    {0x55555555, 0x0763638e},
    {0x2aaaaaaa, 0x1b6d4949},
    {0x55555555, 0x638ef1c7},
    {0x00000000, 0x47638e39},
    {0x2aaaaaaa, 0x5236b6db},
    {0x55555555, 0x4e381c1c},
    {0x7fffffff, 0x638eb1c7},
    {0x7fffffff, 0x47638e39},
    {0x3f2bb31e, 0xb6c004cc}
};

const pair<uint32_t, uint32_t> consts_Malzette2[] = {
    {0x1c71c924, 0x249cad47},
    {0x49249c71, 0x1249871c},
    {0x6db6c71c, 0x5b127ffe},
    {0x38e39249, 0x152ad249},
    {0x638e36db, 0x649cad55},
    {0x1c71c7ff, 0x471c9492},
    {0x36db6d55, 0x63f1c71d},
    {0x471c7249, 0x36a4ff1c},
    {0x4924938e, 0x5b6c8e47},
    {0x2aab6db6, 0x71c736db},
    {0x6db638e3, 0x55b9c71d},
    {0xfb3d2330, 0xb6da4b61}
};

const pair<uint32_t, uint32_t> differential_Malzette1[] = {
    {0x7ffffff8, 0x3ffffffe},
    {0xb989d417, 0x0347a2a9}
};

const pair<uint32_t, uint32_t> differential_Malzette2[] = {
    {0x7fff0003, 0xffffc001},
    {0xdd2bc54b, 0x078b106c}
};

const uint32_t differential_trail_Malzette1[] = {
    0xfffffffc,
    0x7ffffff8,
    0x80000007,
    0x80000007,
    0xfffffffc,
    0x00000003,
    0x7ffffff8,
    0xfffffffc,
    0x00000003,
    0x7ffffff8,
    0x7ffffff8,
    0x7ffffff8,
};
const uint32_t differential_trail_Malzette2[] = {
    0x7fff8007,
    0x80000007,
    0x0000fff8,
    0x80000003,
    0xffff8007,
    0x80000003,
    0xfffff804,
    0xfffffffc,
    0x80000ffb,
    0x00000003,
    0x7ff80004,
    0xfffffffc,
};

const pair<int, int> rotations[] = {
    {31, 24},
    {17, 17},
    {0, 31},
    {24, 16},

    {31, 24},
    {17, 17},
    {0, 31},
    {24, 16},

    {31, 24},
    {17, 17},
    {0, 31},
    {24, 16}
};


static void Malzette(uint32_t &l, uint32_t &r, const pair<uint32_t, uint32_t> *consts) {
	for(int i = 0; i < ROUNDS; i++) {
		l += ror(r, rotations[i].first);
		r ^= ror(l, rotations[i].second);
		l ^= consts[i].first;
		r ^= consts[i].second;
	}
}

static bool followsTrail(uint32_t l, uint32_t r, uint32_t ll, uint32_t rr, const uint32_t *trail, const pair<uint32_t, uint32_t> *consts) {
	for(int i = 0; i < ROUNDS; i++) {
		l += ror(r, rotations[i].first);
		ll += ror(rr, rotations[i].first);

		uint32_t delta = rol(l, K) ^ ll;
		if (delta != trail[i]) {
			printf("  failed follow at round %d/%d: %08x != %08x\n", i+1, ROUNDS, delta, trail[i]);
			return false;
		}

		r ^= ror(l, rotations[i].second);
		rr ^= ror(ll, rotations[i].second);
		
		l ^= consts[i].first;
		ll ^= consts[i].first;
		
		r ^= consts[i].second;
		rr ^= consts[i].second;
	}
	return true;
}

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

int main(int argc, char *argv[]) {
	srandom(time(NULL));
	srandom(random() ^ getpid());
	uint64_t seed1 = (random() << 32) ^ random();
	uint64_t seed2 = (random() << 32) ^ random();

	// 2^34 ~= 2.5 mins on 8 threads
	int data_log = 32;

	if (argc > 1) {
		data_log = atoi(argv[1]);
	}
	if (argc > 3) {
		seed1 = strtoull(argv[2], NULL, 16);
		seed2 = strtoull(argv[3], NULL, 16);
	}

	for(int version = 1; version <= 2; version++) {
		printf("Testing Malzette%d, using 2^%d data, seeds %016lx %016lx\n", version, data_log, seed1, seed2);

		auto differential_Malzette = (version == 1) ? \
			differential_Malzette1 : differential_Malzette2; 
		auto trail_Malzette = (version == 1) ? \
			differential_trail_Malzette1 : differential_trail_Malzette2; 
		auto consts_Malzette = (version == 1) ? \
			consts_Malzette1 : consts_Malzette2;
		
		uint64_t hits = 0;
		uint64_t hits_trail = 0;
		uint64_t n_iters = 1ull << data_log;
		
		#pragma omp parallel for
		for (uint64_t itr = 0; itr < n_iters; itr++) {
			uint32_t l0, r0, l, r, ll, rr, ll0, rr0;
			// note: random() slows openmp a lot!
			// l0 = random() ^ (random() << 16);
			// r0 = random() ^ (random() << 16);;

			// use custom random if openmp
			l0 = hash32(seed1 ^ itr);
			r0 = hash32(seed2 ^ itr);

			l = l0;
			r = r0;
			Malzette(l, r, consts_Malzette);

			ll = ll0 = rol(l0, K) ^ differential_Malzette[0].first;
			rr = rr0 = rol(r0, K) ^ differential_Malzette[0].second;
			Malzette(ll, rr, consts_Malzette);

			if (ll != (rol(l, K) ^ differential_Malzette[1].first))
				continue;
			if (rr != (rol(r, K) ^ differential_Malzette[1].second))
				continue;
			
			#pragma omp critical
			{
				hits++;
				bool follows = followsTrail(l0, r0, ll0, rr0, trail_Malzette, consts_Malzette);
				hits_trail += follows;
				printf("Malzette%d: found match #%lu %08x:%08x (follows trail? %d)\n",
					version, hits, l0, r0, follows
				);
				fflush(stdout);
			}
		};

		printf("Malzette%d: %lu/%lu = 2^%.2f experimental differential probability\n",
			version, hits, n_iters, (log(hits)-log(n_iters))/log(2));
		printf("Malzette%d: %lu/%lu/%lu = 2^%.2f experimental trail probability\n",
			version, hits_trail, hits, n_iters, (log(hits_trail)-log(n_iters))/log(2));
		printf("\n");
		fflush(stdout);
	}
	return 0;
}
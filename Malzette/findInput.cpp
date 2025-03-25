#include <bits/stdc++.h>

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


void Malzette(uint32_t &l, uint32_t &r, const pair<uint32_t, uint32_t> *consts) {
	for(int i = 0; i < ROUNDS; i++) {
		l += ror(r, rotations[i].first);
		r ^= ror(l, rotations[i].second);
		l ^= consts[i].first;
		r ^= consts[i].second;
	}
}


int main() {
	srand(2024);
	uint32_t l0, l, ll, r0, r, rr;
	
	uint64_t itr;
	
	itr = 0;
	while (++itr) {
		l0 = rand();
		r0 = rand();
		if (itr % 16000000 == 0) {
			printf("itr %lu %08x %08x\n", itr, l0, r0);
		}
		
		l = l0;
		r = r0;
		Malzette(l, r, consts_Malzette1);

		ll = rol(l0, K) ^ differential_Malzette1[0].first;
		rr = rol(r0, K) ^ differential_Malzette1[0].second;
		Malzette(ll, rr, consts_Malzette1);

		if (ll != (rol(l, K) ^ differential_Malzette1[1].first))
			continue;
		if (rr != (rol(r, K) ^ differential_Malzette1[1].second))
			continue;
		printf("Malzette1: found match %08x:%08x at iteration %lu = 2^%.2f\n", l0, r0, itr, log(itr)/log(2));
		break;
	};

	itr = 0;
	while (++itr) {
		l0 = rand();
		r0 = rand();
		if (itr % 16000000 == 0) {
			printf("itr %lu %08x %08x\n", itr, l0, r0);
		}
		
		l = l0;
		r = r0;
		Malzette(l, r, consts_Malzette2);

		ll = rol(l0, K) ^ differential_Malzette2[0].first;
		rr = rol(r0, K) ^ differential_Malzette2[0].second;
		Malzette(ll, rr, consts_Malzette2);

		if (ll != (rol(l, K) ^ differential_Malzette2[1].first))
			continue;
		if (rr != (rol(r, K) ^ differential_Malzette2[1].second))
			continue;
		printf("Malzette2: found match %08x:%08x at iteration %lu = 2^%.2f\n", l0, r0, itr, log(itr)/log(2));
		break;
	};
	return 0;
}
#include "common.hpp"

using namespace std;

uint64_t CSHL(uint64_t const x,
		  uint const s,
		  uint const n){
	/*
	Return the cyclic-shift of x by s position to the left, with x over n significant bits
	*/
	if(s == 0) return (x & ((1ULL << n) - 1));
	return ((x << s) & ((1ULL << n) - 1)) | (x >> (n - s));
}

uint64_t CSHR(uint64_t const x,
		  uint const s,
		  uint const n){
	/*
	Return the cyclic-shift of x by s position to the right, with x over n bits
	*/
	if(s == 0) return (x & ((1ULL << n) - 1));
	return (x >> s) | ((x & ((1ULL << s)-1)) << (n-s));
}

uint64_t shl(uint64_t const x,
			 uint const s,
			 uint const n){
	return ((x << s) & ((1ULL << n) - 1));
}

uint64_t shr(uint64_t const x,
			 uint const s,
			 uint const n){
	return ((x >> s) & ((1ULL << n) - 1));
}

void binprint(uint64_t const x,
			  uint const n,
			  bool const rightLSB){
	if(rightLSB){
		for(int i = n-1; i >= 0; i--){
			if(x & (1ULL << i)) cout << "1";
			else cout << "0";
		}
	}
	else{
		for(uint i = 0; i < n; i++){
			if(x & (1ULL << i)) cout << "1";
			else cout << "0";
		}
	}
}

int mod(int const x,
		int const n){
	return ((x%n) + n)%n;
}

uint64_t binVectorToUint(std::vector<uint> const & v){
	uint64_t x = 0;
	for(uint i = 0; i < v.size(); i++){
		if(v[i] != 0)
			x |= (1ULL << i);
	}
	return x;
}

uint64_t nextUintSameHammingWeight(uint const val){
	uint64_t c = val & (-val);
	uint64_t r = val + c;
	return (((r^val) >> 2) / c) | r;
}

std::set<uint64_t> 
genAllVectorHammingDistance(std::vector<std::vector<uint>> const & knownVectors,
							uint const distance){

	set<uint64_t> allVectors;

	for(auto const & v : knownVectors){
		uint64_t base = binVectorToUint(v);
		allVectors.emplace(base);
		uint64_t bound = (1ULL << v.size());

		//For each hamming distance from 1 to max distance
		for(uint d = 1; d <= distance; d++){

			uint64_t mask = (1ULL << d)-1;
			//mask is going to go through all integers with hamming weight d
			while(mask < bound){
				//XOR the mask to generate a binvector with distance d from base
				allVectors.emplace(base^mask);
				//update the mask to the next uint with same hamming weight
				mask = nextUintSameHammingWeight(mask);
			}
		}	
	}
	return allVectors;
}

void writeVectorUint64ToFile(std::vector<uint64_t> const & v,
							 std::string const & filename,
							 std::ios_base::openmode mode){

	ofstream fileout(filename,mode);
	size_t size = v.size();
	fileout.write((char*)&size, sizeof(size));
	fileout.write((char*)&v[0], v.size() * sizeof(uint64_t));
	fileout.close();
}

std::vector<uint64_t> readVectorUint64FromFile(std::string const & filename,
											   std::ios_base::openmode mode,
											   uint64_t const offset){

	ifstream filein(filename,mode);
	filein.seekg(offset);
	size_t size;
	filein.read((char*)&size, sizeof(size));
	vector<uint64_t> v(size);
	filein.read((char*)&v[0], v.size() * sizeof(uint64_t));
	return v;
}

bool checkCriteria(uint64_t const a,
				   uint64_t const b,
				   uint64_t const d,
				   uint const n,
				   uint const k){

	uint64_t mask = (1ULL << n) - 1;
	uint64_t u = (a^b^d) ^ shl(a^b^d,1,n);
	uint64_t v = shl( (a^d) | (b^d) ,1,n);

	if(k == 0){
		return ((u&v) == u);
	}
	else{
		uint64_t mask1 = (~(1ULL))&mask;
		uint64_t maskgamma = (~(1ULL << k))&mask;
		u = u&mask1&maskgamma;
		v = v&mask1&maskgamma;
		return ((u&v) == u);
	}
}


uint64_t Tcounter(uint64_t const alpha,
				  uint64_t const beta,
				  uint64_t const delta,
				  uint const n,
				  uint const v){

	// cout << "Tcounter " << alpha << " " << beta << " " << delta << " " << n << " " << v << endl;
	uint64_t mask = (1ULL << n) - 1;
	uint64_t tmp = (((alpha ^ delta) | (beta ^ delta)) << 1)&mask;
	uint d = __builtin_popcountll(tmp);
	uint64_t xorall = alpha^beta^delta;

	// cout << "tmp " << tmp << endl;
	// cout << "d " << d << endl;
	// cout << "xorall " << xorall << endl;

	if(xorall == 0 || xorall == mask){
		if(v == 0){ 
			return (1ULL << (2*n-d-1)) + (1ULL << (n-1));
		}
		else{
			return (1ULL << (2*n-d-1)) - (1ULL << (n-1));
		}
	}
	else{
		//Number of solutions = 2^(2n-d-1)
		return (1ULL << (2*n-d-1));
	}
}

uint64_t getRXDiffCount(uint64_t const alpha,
						uint64_t const beta,
						uint64_t const delta,
						uint const n,
						uint const k){

	if(!checkCriteria(alpha,beta,delta,n,k))
		return 0;

	if(k == 0){
		uint64_t mask = (1ULL << n) - 1;
		uint64_t tmp = (((alpha ^ delta) | (beta ^ delta)) << 1)&mask;
		uint d = __builtin_popcountll(tmp);
		return (1ULL << (2*n-d));
	}
	else{

		uint64_t mask_k = (1ULL << k)-1;
		uint64_t mask_nk = (1ULL << (n-k))-1;

		uint64_t alpha_L = (alpha >> k)&mask_nk;
		uint64_t  beta_L = (beta  >> k)&mask_nk;
		uint64_t delta_L = (delta >> k)&mask_nk;

		uint64_t alpha_R = alpha&mask_k;
		uint64_t  beta_R = beta&mask_k;
		uint64_t delta_R = delta&mask_k;

		uint v  = (alpha_L^beta_L^delta_L)&1;
		uint vp = (alpha_R^beta_R^delta_R)&1;

		uint64_t factornk = Tcounter(alpha_L,beta_L,delta_L,n-k,vp);
		uint64_t factork = Tcounter(alpha_R,beta_R,delta_R,k,v);

		// cout << "factork  : " << factork << endl;
		// cout << "factornk : " << factornk << endl;

		// uint64_t factornk = 0;
		// if(v == 0)
		// 	factornk = Tcounter(alpha_L,beta_L,delta_L,n-k,vp);
		// else
		// 	factornk = Tcounter(alpha_L^mask_nk,beta_L^mask_nk,delta_L^mask_nk,n-k,vp);

		// uint64_t factork = 0;
		// if(vp == 0)
		// 	factork = Tcounter(alpha_R,beta_R,delta_R,k,v);
		// else
		// 	factork = Tcounter(alpha_R^mask_k,beta_R^mask_k,delta_R^mask_k,k,v);

		return factornk*factork;
	}
}

void checkCountFormula(uint const n, uint const k){

	uint64_t bound_n = (1ULL << n);
	uint64_t mask_n = bound_n - 1;
	vector<vector<vector<uint64_t>>> ddt(bound_n, vector<vector<uint64_t>>(bound_n, vector<uint64_t>(bound_n, 0)));
	for(uint64_t alpha = 0; alpha < bound_n; alpha++){
		for(uint64_t beta = 0; beta < bound_n; beta++){

			for(uint64_t x = 0; x < bound_n; x++){
				for(uint64_t y = 0; y < bound_n; y++){

					uint64_t delta = CSHL((x+y)&mask_n,k,n) ^
									 (((CSHL(x,k,n)^alpha) + (CSHL(y,k,n)^beta))&mask_n);

					ddt[alpha][beta][delta]++;
				}
			}
		}
	}

	uint64_t maxCount = 0;
	uint64_t minCount = -1;
	for(uint64_t alpha = 0; alpha < bound_n; alpha++){
		for(uint64_t beta = 0; beta < bound_n; beta++){
			for(uint64_t delta = 0; delta < bound_n; delta++){
				uint64_t count = getRXDiffCount(alpha,beta,delta,n,k);
				if(ddt[alpha][beta][delta] != count){

					cout << "alpha = "; binprint(alpha,n);
					cout << " beta = "; binprint(beta,n);
					cout << " delta = "; binprint(delta,n);
					cout << " experimental : " << ddt[alpha][beta][delta];
					cout << " theoretical : " << count << endl;
				}
				if(count > maxCount)
					maxCount = count;
				if(count < minCount && count > 0)
					minCount = count;
			}
		}
	}
	cout << "minCount = " << minCount << " maxCount = " << maxCount << endl;
}
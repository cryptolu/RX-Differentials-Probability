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
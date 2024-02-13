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
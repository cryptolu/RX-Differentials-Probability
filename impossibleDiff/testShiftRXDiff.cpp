#include "testShiftRXDiff.hpp"

using namespace std;

uint getBit(uint64_t const x,
			int const i,
			uint const n){
	return (x >> mod(i,n))&1;
}

void checkValueSHL(int const n,
				   int const gamma,
				   int const alpha,
				   uint64_t const x0,
				   uint64_t const x1,
				   uint64_t const deltaOut){

	if(gamma+alpha <= n){
		if(gamma <= alpha){
			stringstream ss;
			ss << "Error: (alpha,gamma) = (" << alpha << "," << gamma << "), case gamma+alpha <= n && gamma <= alpha; i = ";
			string errorPrefix = ss.str();

			for(int i = 0; i < gamma; i++){
				if(getBit(deltaOut,i,n) != getBit(x0,i-gamma-alpha,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma; i < alpha; i++){
				if(getBit(deltaOut,i,n) != 0)
					cout << errorPrefix << i << endl;
			}
			for(int i = alpha; i < gamma+alpha; i++){
				if(getBit(deltaOut,i,n) != getBit(x1,i-alpha,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma+alpha; i < n; i++){
				if(getBit(deltaOut,i,n) != (getBit(x0,i-gamma-alpha,n)^getBit(x1,i-alpha,n)))
					cout << errorPrefix << i << endl;
			}
		}
		else{ //gamma > alpha
			stringstream ss;
			ss << "Error: (alpha,gamma) = (" << alpha << "," << gamma << "), case gamma+alpha <= n && gamma > alpha; i = ";
			string errorPrefix = ss.str();

			for(int i = 0; i < alpha; i++){
				if(getBit(deltaOut,i,n) != getBit(x0,i-gamma-alpha,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = alpha; i < gamma; i++){
				if(getBit(deltaOut,i,n) != (getBit(x0,i-gamma-alpha,n)^getBit(x1,i-alpha,n)))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma; i < gamma+alpha; i++){
				if(getBit(deltaOut,i,n) != getBit(x1,i-alpha,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma+alpha; i < n; i++){
				if(getBit(deltaOut,i,n) != (getBit(x0,i-gamma-alpha,n)^getBit(x1,i-alpha,n)))
					cout << errorPrefix << i << endl;
			}
		}
	}
	else{ //gamma+alpha > n
		if(gamma <= alpha){
			stringstream ss;
			ss << "Error: (alpha,gamma) = (" << alpha << "," << gamma << "), case gamma+alpha > n && gamma <= alpha; i = ";
			string errorPrefix = ss.str();

			for(int i = 0; i < gamma+alpha-n; i++){
				if(getBit(deltaOut,i,n) != 0)
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma+alpha-n; i < gamma; i++){
				if(getBit(deltaOut,i,n) != getBit(x0,i-gamma-alpha,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma; i < alpha; i++){
				if(getBit(deltaOut,i,n) != 0)
					cout << errorPrefix << i << endl;
			}
			for(int i = alpha; i < n; i++){
				if(getBit(deltaOut,i,n) != getBit(x1,i-alpha,n))
					cout << errorPrefix << i << endl;
			}
		}
		else{ //gamma > alpha
			stringstream ss;
			ss << "Error: (alpha,gamma) = (" << alpha << "," << gamma << "), case gamma+alpha > n && gamma > alpha; i = ";
			string errorPrefix = ss.str();

			for(int i = 0; i < gamma+alpha-n; i++){
				if(getBit(deltaOut,i,n) != 0)
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma+alpha-n; i < alpha; i++){
				if(getBit(deltaOut,i,n) != getBit(x0,i-gamma-alpha,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = alpha; i < gamma; i++){
				if(getBit(deltaOut,i,n) != (getBit(x0,i-gamma-alpha,n)^getBit(x1,i-alpha,n)))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma; i < n; i++){
				if(getBit(deltaOut,i,n) != getBit(x1,i-alpha,n))
					cout << errorPrefix << i << endl;
			}
		}
	}
}

void testValueSHL(uint const n){

	cout << "Test value SHL n = " << n << endl;

	uint64_t bound = (1ULL << n);

	for(uint gamma = 1; gamma < n; gamma++){
		for(uint alpha = 1; alpha < n; alpha++){
			cout << "--- (alpha,gamma) = (" << alpha << "," << gamma << ") ---" << endl;
			for(uint64_t x0 = 0; x0 < bound; x0++){
				uint64_t ry0 = CSHL(shl(x0,alpha,n), gamma,n);
				for(uint64_t x1 = 0; x1 < bound; x1++){
					uint64_t y1 = shl(x1,alpha,n);
					uint64_t deltaOut = ry0 ^ y1;
					checkValueSHL(n,gamma,alpha,x0,x1,deltaOut);
				}
			}
		}
	}
}

void checkValueSHR(int const n,
				   int const gamma,
				   int const beta,
				   uint64_t const x0,
				   uint64_t const x1,
				   uint64_t const deltaOut){

	if(gamma+beta <= n){
		if(gamma <= beta){
			stringstream ss;
			ss << "Error: (beta,gamma) = (" << beta << "," << gamma << "), case gamma+beta <= n && gamma <= beta; i = ";
			string errorPrefix = ss.str();

			for(int i = 0; i < gamma; i++){
				if(getBit(deltaOut,i,n) != getBit(x1,i+beta,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma; i < n-beta; i++){
				if(getBit(deltaOut,i,n) != (getBit(x0,i+beta-gamma,n)^getBit(x1,i+beta,n)))
					cout << errorPrefix << i << endl;
			}
			for(int i = n-beta; i < n-beta+gamma; i++){
				if(getBit(deltaOut,i,n) != getBit(x0,i-gamma+beta,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = n-beta+gamma; i < n; i++){
				if(getBit(deltaOut,i,n) != 0)
					cout << errorPrefix << i << endl;
			}
		}
		else{ //gamma > beta
			stringstream ss;
			ss << "Error: (beta,gamma) = (" << beta << "," << gamma << "), case gamma+beta <= n && gamma > beta; i = ";
			string errorPrefix = ss.str();

			for(int i = 0; i < gamma-beta; i++){
				if(getBit(deltaOut,i,n) != (getBit(x0,i+beta-gamma,n)^getBit(x1,i+beta,n)))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma-beta; i < gamma; i++){
				if(getBit(deltaOut,i,n) != getBit(x1,i+beta,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma; i < n-beta; i++){
				if(getBit(deltaOut,i,n) != (getBit(x0,i+beta-gamma,n)^getBit(x1,i+beta,n)))
					cout << errorPrefix << i << endl;
			}
			for(int i = n-beta; i < n; i++){
				if(getBit(deltaOut,i,n) != getBit(x0,i+beta-gamma,n))
					cout << errorPrefix << i << endl;
			}
		}
	}
	else{ //gamma+beta > n
		if(gamma <= beta){
			stringstream ss;
			ss << "Error: (beta,gamma) = (" << beta << "," << gamma << "), case gamma+beta > n && gamma <= beta; i = ";
			string errorPrefix = ss.str();

			for(int i = 0; i < n-beta; i++){
				if(getBit(deltaOut,i,n) != getBit(x1,i+beta,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = n-beta; i < gamma; i++){
				if(getBit(deltaOut,i,n) != 0)
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma; i < gamma+n-beta; i++){
				if(getBit(deltaOut,i,n) != getBit(x0,i+beta-gamma,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma+n-beta; i < n; i++){
				if(getBit(deltaOut,i,n) != 0)
					cout << errorPrefix << i << endl;
			}
		}
		else{ //gamma > beta
			stringstream ss;
			ss << "Error: (beta,gamma) = (" << beta << "," << gamma << "), case gamma+beta > n && gamma > beta; i = ";
			string errorPrefix = ss.str();

			for(int i = 0; i < gamma-beta; i++){
				if(getBit(deltaOut,i,n) != (getBit(x0,i+beta-gamma,n)^getBit(x1,i+beta,n)))
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma-beta; i < n-beta; i++){
				if(getBit(deltaOut,i,n) != getBit(x1,i+beta,n))
					cout << errorPrefix << i << endl;
			}
			for(int i = n-beta; i < gamma; i++){
				if(getBit(deltaOut,i,n) != 0)
					cout << errorPrefix << i << endl;
			}
			for(int i = gamma; i < n; i++){
				if(getBit(deltaOut,i,n) != getBit(x0,i+beta-gamma,n))
					cout << errorPrefix << i << endl;
			}
		}
	}
}

void testValueSHR(uint const n){

	cout << "Test value SHR n = " << n << endl;
	uint64_t bound = (1ULL << n);

	for(uint gamma = 1; gamma < n; gamma++){
		for(uint beta = 1; beta < n; beta++){
			cout << "--- (beta,gamma) = (" << beta << "," << gamma << ") ---" << endl;
			for(uint64_t x0 = 0; x0 < bound; x0++){
				uint64_t ry0 = CSHL(shr(x0,beta,n), gamma,n);
				for(uint64_t x1 = 0; x1 < bound; x1++){
					uint64_t y1 = shr(x1,beta,n);
					uint64_t deltaOut = ry0 ^ y1;
					checkValueSHR(n,gamma,beta,x0,x1,deltaOut);
				}
			}
		}
	}
}

bool checkRXDiffSHL(int const n,
					int const gamma,
					int const alpha,
					uint64_t const deltaIn,
					uint64_t const deltaOut){

	if(gamma+alpha <= n){
		if(gamma <= alpha){
			uint64_t mask_nag = (1ULL << (n-alpha-gamma))-1;
			uint64_t mask_ag = (1ULL << (alpha-gamma))-1;
			if(    ((deltaOut >> gamma)&mask_ag) == 0 
				&& ((deltaOut >> (alpha+gamma))&mask_nag) == ((deltaIn >> gamma)&mask_nag))
				return true;
			else
				return false;
		}
		else{ //gamma > alpha
			uint64_t mask_nag = (1ULL << (n-alpha-gamma))-1;
			uint64_t mask_ga = (1ULL << (gamma-alpha))-1;
			if(    ((deltaOut >> alpha)&mask_ga) == (deltaIn&mask_ga) 
				&& ((deltaOut >> (alpha+gamma))&mask_nag) == ((deltaIn >> gamma)&mask_nag))
				return true;
			else
				return false;
		}
	}
	else{ //gamma+alpha > n
		if(gamma <= alpha){
			uint64_t mask_gan = (1ULL << (gamma+alpha-n))-1;
			uint64_t mask_ag = (1ULL << (alpha-gamma))-1;
			if(    (deltaOut&mask_gan) == 0
				&& ((deltaOut >> gamma)&mask_ag) == 0)
				return true;
			else
				return false;
		}
		else{ //gamma > alpha
			uint64_t mask_gan = (1ULL << (gamma+alpha-n))-1;
			uint64_t mask_ga = (1ULL << (gamma-alpha))-1;
			if(    (deltaOut&mask_gan) == 0
				&& ((deltaOut >> alpha)&mask_ga) == (deltaIn&mask_ga))
				return true;
			else
				return false;
		}
	}
}

void testRXDiffSHL(uint const n){

	cout << "Test RX-diff SHL n = " << n << endl;

	uint64_t bound = (1ULL << n);

	for(uint gamma = 1; gamma < n; gamma++){
		for(uint alpha = 1; alpha < n; alpha++){
			cout << "--- (alpha,gamma) = (" << alpha << "," << gamma << ") ---" << endl;

			//precompute the full table
			vector<vector<bool>> table(bound, vector<bool>(bound,false));
			for(uint64_t x0 = 0; x0 < bound; x0++){
				uint64_t rx0 = CSHL(x0, gamma,n);
				uint64_t ry0 = CSHL(shl(x0,alpha,n), gamma,n);
				for(uint64_t x1 = 0; x1 < bound; x1++){
					uint64_t y1 = shl(x1,alpha,n);
					uint64_t deltaIn = rx0 ^ x1;
					uint64_t deltaOut = ry0 ^ y1;
					table[deltaIn][deltaOut] = true;
				}
			}

			//Check that everything is consistent
			for(uint64_t deltaIn = 0; deltaIn < bound; deltaIn++){
				for(uint64_t deltaOut = 0; deltaOut < bound; deltaOut++){
					if(checkRXDiffSHL(n,gamma,alpha,deltaIn,deltaOut) != table[deltaIn][deltaOut]){
						cout << "Error on deltaIn = " << deltaIn << " deltaOut = " << deltaOut << endl;
						exit(1);
					}
				}
			}
		}
	}
}

bool checkRXDiffSHR(int const n,
					int const gamma,
					int const beta,
					uint64_t const deltaIn,
					uint64_t const deltaOut){
	if(gamma+beta <= n){
		if(gamma <= beta){
			uint64_t mask_nbg = (1ULL << (n-beta-gamma))-1;
			uint64_t mask_bg = (1ULL << (beta-gamma))-1;
			if(    ((deltaOut >> gamma)&mask_nbg) == ((deltaIn >> (beta+gamma))&mask_nbg) 
				&& ((deltaOut >> (n+gamma-beta))&mask_bg) == 0)
				return true;
			else
				return false;
		}
		else{ //gamma > beta
			uint64_t mask_nbg = (1ULL << (n-beta-gamma))-1;
			uint64_t mask_gb = (1ULL << (gamma-beta))-1;
			if(    (deltaOut&mask_gb) == ((deltaIn >> beta)&mask_gb) 
				&& ((deltaOut >> gamma)&mask_nbg) == ((deltaIn >> (beta+gamma))&mask_nbg))
				return true;
			else
				return false;
		}
	}
	else{ //gamma+beta > n
		if(gamma <= beta){
			uint64_t mask_gbn = (1ULL << (gamma+beta-n))-1;
			uint64_t mask_bg = (1ULL << (beta-gamma))-1;
			if(    ((deltaOut >> (n-beta))&mask_gbn) == 0
				&& ((deltaOut >> (n+gamma-beta))&mask_bg) == 0)
				return true;
			else
				return false;
		}
		else{ //gamma > beta
			uint64_t mask_gb = (1ULL << (gamma-beta))-1;
			uint64_t mask_gbn = (1ULL << (gamma+beta-n))-1;
			if(    (deltaOut&mask_gb) == ((deltaIn >> beta)&mask_gb) 
				&& ((deltaOut >> (n-beta))&mask_gbn) == 0)
				return true;
			else
				return false;
		}
	}
}

void testRXDiffSHR(uint const n){

	cout << "Test RX-diff SHR n = " << n << endl;

	uint64_t bound = (1ULL << n);

	for(uint gamma = 1; gamma < n; gamma++){
		for(uint beta = 1; beta < n; beta++){
			cout << "--- (beta,gamma) = (" << beta << "," << gamma << ") ---" << endl;

			//precompute the full table
			vector<vector<bool>> table(bound, vector<bool>(bound,false));
			for(uint64_t x0 = 0; x0 < bound; x0++){
				uint64_t rx0 = CSHL(x0, gamma,n);
				uint64_t ry0 = CSHL(shr(x0,beta,n), gamma,n);
				for(uint64_t x1 = 0; x1 < bound; x1++){
					uint64_t y1 = shr(x1,beta,n);
					uint64_t deltaIn = rx0 ^ x1;
					uint64_t deltaOut = ry0 ^ y1;
					table[deltaIn][deltaOut] = true;
				}
			}

			//Check that everything is consistent
			for(uint64_t deltaIn = 0; deltaIn < bound; deltaIn++){
				for(uint64_t deltaOut = 0; deltaOut < bound; deltaOut++){
					if(checkRXDiffSHR(n,gamma,beta,deltaIn,deltaOut) != table[deltaIn][deltaOut]){
						cout << "Error on deltaIn = " << deltaIn << " deltaOut = " << deltaOut << endl;
						exit(1);
					}
				}
			}
		}
	}
}

int main(){
	testValueSHL(8);
	testValueSHR(8);

	testRXDiffSHL(8);
	testRXDiffSHR(8);
}
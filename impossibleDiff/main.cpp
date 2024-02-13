#include <vector>
#include <utility>
#include <cstdint>

#include <omp.h>

#include "speck.hpp"
#include "alzette.hpp"
#include "xtea.hpp"
#include "lea.hpp"
#include "common.hpp"

using namespace std;



int main(){

	//Most of the searches are implemented with some multithreading
	//In most cases, it is a lot of fast to solve models, so dedicating only a small number of threads to Gurobi and processing more models at once through multithreading is more beneficial
	omp_set_num_threads(8);

	/* -----    Search for Alzette truncated RX-differentials   ----- */
	/* ----- This is used to provide the results in Section 6.4 ----- */
	{
	vector<uint> paramsAlzette({31,24,17,17,0,31,24,16});
	uint32_t cstAlzette = 0x4f7c7b57;
	uint gamma = 7;
	vector<uint> inputIndex({34,35});
	vector<uint> outputIndex({12,21,37,60});
	checkTruncatedDiffAlzette(paramsAlzette,cstAlzette,gamma,inputIndex,outputIndex);
	cout << endl << endl;
	}

	// cstAlzette = 0xbf715880;
	// gamma = 1;
	// inputIndex = vector<uint>({63});
	// outputIndex = vector<uint>({12,60});
	// checkTruncatedDiffAlzette(paramsAlzette,cstAlzette,gamma,inputIndex,outputIndex);
	// cout << endl << endl;

	// cstAlzette = 0x324e7738;
	// gamma = 5;
	// inputIndex = vector<uint>({4});
	// outputIndex = vector<uint>({23,39});
	// checkTruncatedDiffAlzette(paramsAlzette,cstAlzette,gamma,inputIndex,outputIndex);
	// cout << endl << endl;

	// cstAlzette = 0xcfbfa1c8;
	// gamma = 3;
	// inputIndex = vector<uint>({31,62});
	// outputIndex = vector<uint>({34});
	// checkTruncatedDiffAlzette(paramsAlzette,cstAlzette,gamma,inputIndex,outputIndex);
	// cout << endl << endl;

	// cstAlzette = 0xcfbfa1c8;
	// gamma = 3;
	// inputIndex = vector<uint>({34});
	// outputIndex = vector<uint>({12,21,37,60});
	// checkTruncatedDiffAlzette(paramsAlzette,cstAlzette,gamma,inputIndex,outputIndex);

	// cstAlzette = 0xcfbfa1c8;
	// gamma = 29;
	// inputIndex = vector<uint>({63});
	// outputIndex = vector<uint>({9,18,34,57});
	// checkTruncatedDiffAlzette(paramsAlzette,cstAlzette,gamma,inputIndex,outputIndex);

	
	/* ----- Search for Speck related-key truncated RX-differentials   ----- */
	/* -----    This is used to provide the results in Section 6.3     ----- */

	uint gamma = 1;
	vector<uint> inputIndex({32,33,35,36});
	vector<uint> outputIndex({15,16,19,20});
	vector<uint> keyIndex({39});
	checkTruncatedDiffSpeckRelatedKey(32,2,3,gamma,inputIndex,outputIndex,keyIndex);


	/* ----- Search for Alzette weight-(1,1) RX-differentials   ----- */

	// vector<uint> paramsAlzette({31,24,17,17,0,31,24,16});
	// uint32_t cstAlzette = 0xb7e15162;
	// searchWeight1RXDiffAlzette(paramsAlzette,cstAlzette);
	// cstAlzette = 0xbf715880;
	// searchWeight1RXDiffAlzette(paramsAlzette,cstAlzette);
	// cstAlzette = 0x38b4da56;
	// searchWeight1RXDiffAlzette(paramsAlzette,cstAlzette);
	// cstAlzette = 0x324e7738;
	// searchWeight1RXDiffAlzette(paramsAlzette,cstAlzette);
	// cstAlzette = 0xbb1185eb;
	// searchWeight1RXDiffAlzette(paramsAlzette,cstAlzette);
	// uint cstAlzette = 0x4f7c7b57;
	// searchWeight1RXDiffAlzette(paramsAlzette,cstAlzette);
	// cstAlzette = 0xcfbfa1c8;
	// searchWeight1RXDiffAlzette(paramsAlzette,cstAlzette);
	// cstAlzette = 0xc2b3293d;
	// searchWeight1RXDiffAlzette(paramsAlzette,cstAlzette);


	//uint maxRound = 6;
	//vector<uint> cutIndex; //key bit elimination to try to find weak-key classes
	/* ----- Search for Speck related-key weight-(1,1) RX-differentials   ----- */

	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeckRelatedKey(16,2,r);
	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeckRelatedKey(32,2,r);

	/* ----- Search for Single-Key weight-(1,1) RX-differentials   ----- */

	// //------- XTEA without cutIndex --------
	// for(uint r = 5; r < 7; r++)
	// 	searchWeight1RXDiffXtea(r,cutIndex);


	// //------ LEA without cutIndex ------
	// for(uint r = 2; r < maxRound; r++)
	// 	searchWeight1RXDiffLea(128,r,cutIndex);
	
	// for(uint r = 2; r < maxRound; r++)
	// 	searchWeight1RXDiffLea(192,r,cutIndex);

	// for(uint r = 2; r < maxRound; r++)
	// 	searchWeight1RXDiffLea(256,r,cutIndex);
	

	// //------- Speck without cutIndex --------
	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeck(16,2,r,cutIndex);

	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeck(24,2,r,cutIndex);

	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeck(32,2,r,cutIndex);

	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeck(64,2,r,cutIndex);


	/* ----- Search for weak-Kky weight-(1,1) RX-differentials   ----- */
	/* -----    cutIndex set to the first 8 bits of the key      ----- */
	// cutIndex = vector<uint>({0,1,2,3,4,5,6,7,8});

	// //------- XTEA with cutIndex --------
	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffXtea(r,cutIndex);

	// //------ LEA with cutIndex ------
	// for(uint r = 2; r < 5; r++)
	// 	searchWeight1RXDiffLea(128,r,cutIndex);
	
	// for(uint r = 2; r < 5; r++)
	// 	searchWeight1RXDiffLea(192,r,cutIndex);
	
	// for(uint r = 2; r < 5; r++)
	// 	searchWeight1RXDiffLea(256,r,cutIndex);


	// //------- Speck with cutIndex --------
	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeck(16,2,r,cutIndex);

	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeck(24,2,r,cutIndex);

	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeck(32,2,r,cutIndex);

	// for(uint r = 3; r < maxRound; r++)
	// 	searchWeight1RXDiffSpeck(64,2,r,cutIndex);
}
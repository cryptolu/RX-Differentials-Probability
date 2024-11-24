#include <vector>
#include <utility>
#include <cstdint>
#include <random>
#include <fstream>
#include <set>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <assert.h>

#include <omp.h>

#include "ciphers/alzette.hpp"
#include "common/common.hpp"
#include "common/customCallback.hpp"

using namespace std;

const int stdoutfd(dup(fileno(stdout)));

int constIndex = 0;
int rotationAmount = 0;
bool fullOpt = false;
double timeLimitNEQ = 900.0;
double timeLimitFull = 3600.0;

const uint n = 32;
const vector<uint64_t> constAlzette({0xb7e15162,0xbf715880,0x38b4da56,0x324e7738,0xbb1185eb,0x4f7c7b57,0xcfbfa1c8,0xc2b3293d});
const vector<uint> roundRotations({31,24,17,17,0,31,24,16});

void alzetteSearch(){
	GRBEnv env = GRBEnv();
	//----------------------
	//--- Alzette Search ---
	//----------------------
	
	//Search for best trail without input/output constraints

	uint64_t cst = constAlzette[constIndex];
	uint rotationOffset = rotationAmount;
	if (getenv("NO_REDIRECT") == NULL)
	{
		char log_filename[4096];
		sprintf(log_filename,
			"logs/alzette_const%08lx_rotation%u_%s_time%lu.txt",
			cst, rotationOffset,
			(fullOpt ? "fullopt" : "neqopt"),
			(uint64_t)(fullOpt ? timeLimitFull : timeLimitNEQ)
		);
		fflush(stdout);
		fflush(stderr);
		int newstdout = open(log_filename, O_TRUNC | O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
		assert(newstdout > 0);
		dup2(newstdout, fileno(stdout));
		dup2(newstdout, fileno(stderr));
		close(newstdout);
	}

	
	char filename[4096];
	cout << "----------------------------------------------------" << endl;
	cout << "--- Alzette 4 rounds with k = " << rotationOffset;
	cout << " cst = " << "0x" << setfill('0') << setw(n/4) << hex << cst << dec;
	cout << " ";
	cout << (fullOpt ? "full optimization" : "NEQ optimization");
	cout << " TL ";
	cout << timeLimitNEQ << timeLimitFull;
	cout << " ";
	cout << " ---" << endl;
	cout << "----------------------------------------------------" << endl;

	
	
	sprintf(filename,
		"trails/alzette_const%08lx_rotation%u_neqopt_time%lu.txt",
		cst, rotationOffset, (uint64_t)timeLimitNEQ);
	if (!fullOpt) {
		setFilename(filename);

		env.set(GRB_DoubleParam_TimeLimit, timeLimitNEQ);
		searchBestTrailAlzette(roundRotations,cst,rotationOffset,env,false);
		return;
	}
	
	Trail knownTrail;
	knownTrail.read_from_file(filename);
	printf("read trail from file: neq %u log %lf\n", knownTrail.sumNeq, knownTrail.sumLog);
	fflush(stdout);

	if (fullOpt) {
		sprintf(filename,
			"trails/alzette_const%08lx_rotation%u_fullopt_time%lu.txt",
			cst, rotationOffset, (uint64_t)timeLimitFull);
		setFilename(filename);

		env.set(GRB_DoubleParam_TimeLimit, timeLimitFull);
		Trail finalTrail = searchBestTrailAlzette(roundRotations,cst,rotationOffset,env,true,knownTrail);
		
		setFilename(NULL);
		env.set(GRB_DoubleParam_TimeLimit, GRB_INFINITY); // TL = inf
		findPlaintextForTrailAlzette(roundRotations,cst,rotationOffset,env,finalTrail);
	}
}

int main(int argc, char *argv[]){
	if (argc != 1 + 5) {
		printf("Usage: %s <constIndex> <rotationAmount> <fullOpt-flag> <timeLimitNEQ> <timeLimitFull>\nTimelimits in seconds.\n", argv[0]);
		return -1;
	}

	constIndex = atoi(argv[1]);
	rotationAmount = atoi(argv[2]);
	fullOpt = atoi(argv[3]);
	sscanf(argv[4], "%lf", &timeLimitNEQ);
	sscanf(argv[5], "%lf", &timeLimitFull);

	fprintf(stderr, "%s %s %s\n", argv[1], argv[2], argv[3]);
	fprintf(stderr, "%d %d %d\n", constIndex, rotationAmount, fullOpt);
	fflush(stderr);
	assert(0 <= constIndex && (size_t)constIndex < constAlzette.size());
	assert(1 <= rotationAmount && (uint)rotationAmount <= n/2);
	
	alzetteSearch();


/*
	vector<uint64_t> constAlzette({0xb7e15162,0xbf715880,0x38b4da56,0x324e7738,0xbb1185eb,0x4f7c7b57,0xcfbfa1c8,0xc2b3293d});
	
	GRBEnv env = GRBEnv();
	// env.set(GRB_IntParam_OutputFlag, 0);

	uint rotationOffset = 1;
	uint n = 32;

	//4-round trail for Alzette(c[0])
	uint64_t input_x  = 0x0841FA08;
	uint64_t input_y  = 0x0420F900;
	uint64_t output_x = 0xDC33F049;
	uint64_t output_y = 0xBFEE1413;
	vector<uint> params({31,24,17,17,0,31,24,16});
	uint64_t cst = constAlzette[0];

	// //3-round trail for Alzette(c[1])
	// uint64_t input_x  = 0x01000000;
	// uint64_t input_y  = 0x00000000;
	// uint64_t output_x = 0xD5530889;
	// uint64_t output_y = 0x3D412318;
	// vector<uint> params({31,24,17,17,0,31});
	// uint64_t cst = constAlzette[1];

	// vector<uint> inputDiff(2*n);
	vector<uint> outputDiff(2*n);

	// for(uint i = 0; i < n; i++){
	// 	if(((input_x >> i)&1) != 0)
	// 		inputDiff[i] = 1;
	// 	else
	// 		inputDiff[i] = 0;

	// 	if(((input_y >> i)&1) != 0)
	// 		inputDiff[i+n] = 1;
	// 	else
	// 		inputDiff[i+n] = 0;
	// }

	for(uint i = 0; i < n; i++){
		if(((output_x >> i)&1) != 0)
			outputDiff[i] = 1;
		else
			outputDiff[i] = 0;

		if(((output_y >> i)&1) != 0)
			outputDiff[i+n] = 1;
		else
			outputDiff[i+n] = 0;
	}

	vector<uint> inputDiff;
	// vector<uint> outputDiff;
	existTrailAlzette(params,cst,rotationOffset,inputDiff,outputDiff,env);
*/
	return 0;
}
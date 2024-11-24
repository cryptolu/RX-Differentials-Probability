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

#include <omp.h>

#include "alzette.hpp"
#include "speck.hpp"
#include "common.hpp"
#include "salsa.hpp"

using namespace std;



void testTproba(int const n){

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	GRBModel model = GRBModel(env);
	uint64_t bound_n = (1ULL << n);
	double epsilon = 0.000000000000001; //15-digit precision check
	cout << setprecision(16);

	for(uint64_t val_lsb = 0; val_lsb < 2; val_lsb++){
		for(uint64_t val_alpha = 0; val_alpha < bound_n; val_alpha++){
			for(uint64_t val_beta = 0; val_beta < bound_n; val_beta++){
				for(uint64_t val_delta = 0; val_delta < bound_n; val_delta++){

					vector<GRBVar> alpha(n);
					vector<GRBVar> beta(n);
					vector<GRBVar> delta(n);
					for(int i = 0; i < n; i++){
						alpha[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "alpha"+to_string(i));
						beta[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "beta"+to_string(i));
						delta[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "delta"+to_string(i));
					}
					GRBVar lsb = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "lsb");
					GRBVar wt = model.addVar(0,2*n,0,GRB_CONTINUOUS, "wt");
					model.update();

					//differential value constraints
					for(int i = 0; i < n; i++){
						if(((val_alpha >> i)&1) == 0)
							model.addConstr(alpha[i] == 0);
						else
							model.addConstr(alpha[i] == 1);

						if(((val_beta >> i)&1) == 0)
							model.addConstr(beta[i] == 0);
						else
							model.addConstr(beta[i] == 1);

						if(((val_delta >> i)&1) == 0)
							model.addConstr(delta[i] == 0);
						else
							model.addConstr(delta[i] == 1);
					}
					model.addConstr(lsb == val_lsb);

					addTnConstr(model, alpha, beta, delta, lsb, wt, "");

					model.optimize();

					if(model.get(GRB_IntAttr_SolCount) == 0)
						cout << "No solution..." << endl;
					else{
						double expectedlog = -(log2(Tcounter(val_alpha,val_beta,val_delta,n,val_lsb)) - 2*n);
						double resultlog = wt.get(GRB_DoubleAttr_X);

						// if(expectedlog != resultlog){
						if( abs(expectedlog - resultlog) > epsilon){
							cout << val_alpha << " " << val_beta << " " << val_delta << " ";
							cout << val_lsb << " : ";
							cout << resultlog << " ";
							cout << " expected : " << expectedlog << endl;
						}
					}
				}
			}
		}
	}

}

void testModAddProba(int const n, int const k){
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	double epsilon = 0.000000000000001; //15-digit precision check
	cout << setprecision(16);
	cout << "precision of equality check up to epsilon = " << epsilon << endl;
	

	for(uint64_t iter = 0; iter < 1000; iter++){

		std::random_device rd;
	    std::mt19937 gen(rd());

		
		uint64_t val_alpha = 0;
		uint64_t val_beta = 0;
		uint64_t val_delta = 0;
		std::uniform_int_distribution<uint64_t> dis(0, (1ULL << n)-1);
		do{
			val_alpha = dis(gen);
			val_beta = dis(gen);
			val_delta = dis(gen);

		}while(!checkCriteria(val_alpha,val_beta,val_delta,n,k));

		GRBModel model = GRBModel(env);
		//inputs
		vector<GRBVar> alpha(n);
		vector<GRBVar> beta(n);
		vector<GRBVar> delta(n);
		for(int i = 0; i < n; i++){
			alpha[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "alpha"+to_string(i));
			beta[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "beta"+to_string(i));
			delta[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "delta"+to_string(i));
		}
		//output
		GRBVar wt = model.addVar(0,2*n,0,GRB_CONTINUOUS, "wt");
		model.update();

		//differential value constraints
		for(int i = 0; i < n; i++){
			if(((val_alpha >> i)&1) == 0)
				model.addConstr(alpha[i] == 0);
			else
				model.addConstr(alpha[i] == 1);

			if(((val_beta >> i)&1) == 0)
				model.addConstr(beta[i] == 0);
			else
				model.addConstr(beta[i] == 1);

			if(((val_delta >> i)&1) == 0)
				model.addConstr(delta[i] == 0);
			else
				model.addConstr(delta[i] == 1);
		}

		addModAddRXProbaConstr(model, alpha, beta, delta, wt, k, "");

		model.optimize();

		double resultWeight = wt.get(GRB_DoubleAttr_X);
		
		uint64_t expectedCount = getRXDiffCount(val_alpha,val_beta,val_delta,n,k);
		double expectedWeight = -log2(double(expectedCount)/(1ULL << (2*n)));
		
		// if(expectedWeight != resultWeight){
		if( abs(expectedWeight - resultWeight) > epsilon){
			cout << "different weight" << endl;
			cout << "wt           : " << resultWeight << endl;
			cout << "expected log : " << expectedWeight << endl;
		}
	}
}

void alzetteSearch(double timeLimitNeq,
				   double timeLimitFinalTrail,
				   double timeLimitPlaintext){

	GRBEnv env = GRBEnv();
	//----------------------
	//--- Alzette Search ---
	//----------------------
	//Search for best trail without input/output constraints
	vector<uint64_t> cstAlzette({0xb7e15162,0xbf715880,0x38b4da56,0x324e7738,0xbb1185eb,0x4f7c7b57,0xcfbfa1c8,0xc2b3293d});
	vector<uint> params({31,24,17,17,0,31,24,16});
	uint n = 32;

	for(uint cst_index = 0; cst_index < cstAlzette.size(); cst_index++){
		uint64_t cst = cstAlzette[cst_index];
		for(uint rotationOffset = 1; rotationOffset < n/2; rotationOffset++){
			cout << "----------------------------------------------------" << endl;
			cout << "--- Alzette 4 rounds with k = " << rotationOffset;
			cout << " cst = " << "0x" << setfill('0') << setw(n/4) << hex << cst << dec;
			cout << " ---" << endl;
			cout << "----------------------------------------------------" << endl;
			env.set(GRB_DoubleParam_TimeLimit, timeLimitNeq);
			auto knownTrail = searchBestTrailAlzette(params,cst,rotationOffset,env,false);
			env.set(GRB_DoubleParam_TimeLimit, timeLimitFinalTrail);
			auto finalTrail = searchBestTrailAlzette(params,cst,rotationOffset,env,true,knownTrail);
			env.set(GRB_DoubleParam_TimeLimit, timeLimitPlaintext);
			findPlaintextForTrailAlzette(params,cst,rotationOffset,env,finalTrail.second);
		}
	}
}

void speckRelatedKeySearch(uint const n,
						   uint const m,
						   uint const maxRound,
						   double timeLimitNeq,
						   double timeLimitFinalTrail,
						   double timeLimitPlaintext){

	GRBEnv env = GRBEnv();
	// uint n = 16;
	// uint m = 2;
	uint alpha = 8;
	uint beta = 3;
	if(n == 16){
		alpha = 7;
		beta = 2;
	}

	//-------------------------
	//--- Speck Related Key ---
	//-------------------------
	for(uint rotationOffset = 1; rotationOffset < n/2; rotationOffset++){
		map<uint,uint> lowerBoundConsecutiveRounds;
		lowerBoundConsecutiveRounds[1] = 0;
		for(uint nbRound = 2; nbRound <= maxRound; nbRound++){
			array<uint,5> params({n,m,nbRound,alpha,beta});

			cout << "----------------------------------------------------" << endl;
			cout << "--- Speck " << 2*n << "/" << m*n << " " << nbRound << " rounds with k = " << rotationOffset;
			cout << " ---" << endl;
			cout << "----------------------------------------------------" << endl;

			//First heuristic, minimize the number of neq vars
			env.set(GRB_DoubleParam_TimeLimit, timeLimitNeq);
			auto [bestLB, bestNeq, trail] = searchBestTrailSpeckRelatedKey(params,rotationOffset,env,false,lowerBoundConsecutiveRounds);

			lowerBoundConsecutiveRounds[nbRound] = bestLB;
			auto knownTrail = make_pair(bestNeq,trail);

			//Now optimize the probability 
			env.set(GRB_DoubleParam_TimeLimit, timeLimitFinalTrail);
			auto [bestLBproba, bestNeqProba, finalTrail] = searchBestTrailSpeckRelatedKey(params,rotationOffset,env,true,lowerBoundConsecutiveRounds,knownTrail);

			//Find a pair of plaintext to check that the trail is feasible
			env.set(GRB_DoubleParam_TimeLimit, timeLimitPlaintext);
			findPlaintextForTrailSpeckRelatedKey(params,rotationOffset,env,finalTrail);
			getchar();
		}
	}
}

void salsaSearch(uint const maxRound,
				 double timeLimitNeq,
				 double timeLimitFinalTrail,
				 double timeLimitPlaintext){

	// -------------------------
	// ---       Salsa       ---
	// -------------------------
	GRBEnv env = GRBEnv();
	uint startingRound = 0;
	uint n = 32;
	for(uint rotationOffset = 1; rotationOffset < n/2; rotationOffset++){
		map<uint,uint> lowerBoundConsecutiveRounds;
		// lowerBoundConsecutiveRounds[1] = 0;
		for(uint nbRound = 1; nbRound <= maxRound; nbRound++){

			cout << "----------------------------------------------------" << endl;
			cout << "--- Salsa " << nbRound << " rounds with k = " << rotationOffset;
			cout << " ---" << endl;
			cout << "----------------------------------------------------" << endl;

			//First heuristic, minimize the number of neq vars
			env.set(GRB_DoubleParam_TimeLimit, timeLimitNeq);
			env.set(GRB_IntParam_LazyConstraints, 1);
			auto [bestLB, bestNeq, trail] = searchBestTrailSalsa(nbRound,startingRound,rotationOffset,env,false,lowerBoundConsecutiveRounds);

			lowerBoundConsecutiveRounds[nbRound] = bestLB;
			auto knownTrail = make_pair(bestNeq,trail);

			//Now optimize the probability 
			env.set(GRB_DoubleParam_TimeLimit, timeLimitFinalTrail);
			env.set(GRB_IntParam_LazyConstraints, 1);
			auto [bestLBproba, bestNeqProba, finalTrail] = searchBestTrailSalsa(nbRound,startingRound,rotationOffset,env,true,lowerBoundConsecutiveRounds,knownTrail);

			//Find a pair of plaintext to check that the trail is feasible
			env.set(GRB_DoubleParam_TimeLimit, timeLimitPlaintext);
			env.set(GRB_IntParam_LazyConstraints, 0);
			// findPlaintextForTrailSalsa(nbRound,startingRound,rotationOffset,env,finalTrail);

			getchar();
		}
	}

}


int main(){

	double timeLimitNeq = 1800;
	double timeLimitFinalTrail = 1800;
	double timeLimitPlaintext = GRB_INFINITY;

	// //Alzette
	// alzetteSearch(timeLimitNeq,timeLimitFinalTrail,timeLimitPlaintext);

	// //Speck RK
	// uint n = 16;
	// uint m = 2;
	// uint maxRound = 6;
	// speckRelatedKeySearch(n,m,maxRound,timeLimitNeq,timeLimitFinalTrail,timeLimitPlaintext);


	//Salsa
	uint maxRound = 5;
	salsaSearch(maxRound,timeLimitNeq,timeLimitFinalTrail,timeLimitPlaintext);


	// -------------------------
	// ---     Salsa QR      ---
	// -------------------------
	// GRBEnv env = GRBEnv();
	// uint n = 32;
	// for(uint rotationOffset = 1; rotationOffset < n/2; rotationOffset++){
	// 	for(int indexCst = 0; indexCst < 4; indexCst++){

	// 		cout << "----------------------------------------------------" << endl;
	// 		cout << "--- Salsa QR cst " << indexCst << " with k = " << rotationOffset;
	// 		cout << " ---" << endl;
	// 		cout << "----------------------------------------------------" << endl;

	// 		//First heuristic, minimize the number of neq vars
	// 		env.set(GRB_DoubleParam_TimeLimit, timeLimitNeq);
	// 		env.set(GRB_IntParam_LazyConstraints, 1);
	// 		auto [bestLB, bestNeq, trail] = searchBestTrailSalsaQR(indexCst,rotationOffset,env,false);
	// 		auto knownTrail = make_pair(bestNeq,trail);

	// 		//Now optimize the probability 
	// 		env.set(GRB_DoubleParam_TimeLimit, timeLimitFinalTrail);
	// 		env.set(GRB_IntParam_LazyConstraints, 1);
	// 		auto [bestLBproba, bestNeqProba, finalTrail] = searchBestTrailSalsaQR(indexCst,rotationOffset,env,true,knownTrail);

	// 		//Find a pair of plaintext to check that the trail is feasible
	// 		env.set(GRB_DoubleParam_TimeLimit, timeLimitPlaintext);
	// 		env.set(GRB_IntParam_LazyConstraints, 0);
	// 		findPlaintextForTrailSalsaQR(indexCst,rotationOffset,env,finalTrail);
	// 	}
	// }








	

	// //------------------------
	// //--- Speck Single Key ---
	// //------------------------
	// for(uint rotationOffset = 1; rotationOffset < n/2; rotationOffset++){
	// 	map<uint,uint> lowerBoundConsecutiveRounds;
	// 	lowerBoundConsecutiveRounds[1] = 0;
	// 	for(uint nbRound = 2; nbRound < 11; nbRound++){
	// 		array<uint,5> params({n,m,nbRound,alpha,beta});

	// 		cout << "----------------------------------------------------" << endl;
	// 		cout << "--- Speck " << 2*n << "/" << m*n << " " << nbRound << " rounds with k = " << rotationOffset;
	// 		cout << " ---" << endl;
	// 		cout << "----------------------------------------------------" << endl;

	// 		//First heuristic, minimize the number of neq vars
	// 		env.set(GRB_DoubleParam_TimeLimit, timeLimitNeq);
	// 		auto [bestLB, bestNeq, trail] = searchBestTrailSpeck(params,rotationOffset,env,false,lowerBoundConsecutiveRounds);

	// 		lowerBoundConsecutiveRounds[nbRound] = bestLB;
	// 		auto knownTrail = make_pair(bestNeq,trail);

	// 		//Now optimize the probability 
	// 		env.set(GRB_DoubleParam_TimeLimit, timeLimitFinalTrail);
	// 		auto [bestLBproba, bestNeqProba, finalTrail] = searchBestTrailSpeck(params,rotationOffset,env,true,lowerBoundConsecutiveRounds,knownTrail);

	// 		//Find a pair of plaintext to check that the trail is feasible
	// 		env.set(GRB_DoubleParam_TimeLimit, timeLimitPlaintext);
	// 		findPlaintextForTrailSpeck(params,rotationOffset,env,finalTrail);
	// 		getchar();
	// 	}
	// }


	


/*
	vector<uint64_t> cstAlzette({0xb7e15162,0xbf715880,0x38b4da56,0x324e7738,0xbb1185eb,0x4f7c7b57,0xcfbfa1c8,0xc2b3293d});
	
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
	uint64_t cst = cstAlzette[0];

	// //3-round trail for Alzette(c[1])
	// uint64_t input_x  = 0x01000000;
	// uint64_t input_y  = 0x00000000;
	// uint64_t output_x = 0xD5530889;
	// uint64_t output_y = 0x3D412318;
	// vector<uint> params({31,24,17,17,0,31});
	// uint64_t cst = cstAlzette[1];

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
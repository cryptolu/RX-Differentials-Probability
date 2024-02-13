#include "xtea.hpp"
using namespace std;

void addXteaKSConstr(GRBModel & model,
					std::vector<std::vector<GRBVar>> & mk,
					std::vector<std::vector<GRBVar>> & k){

	vector<uint32_t> roundConstants({0x0,0x9e3779b9,0x9e3779b9,0x3c6ef372,0x3c6ef372,0xdaa66d2b,0xdaa66d2b,0x78dde6e4,0x78dde6e4,0x1715609d,0x1715609d,0xb54cda56,0xb54cda56,0x5384540f,0x5384540f,0xf1bbcdc8,0xf1bbcdc8,0x8ff34781,0x8ff34781,0x2e2ac13a,0x2e2ac13a,0xcc623af3,0xcc623af3,0x6a99b4ac,0x6a99b4ac,0x8d12e65,0x8d12e65,0xa708a81e,0xa708a81e,0x454021d7,0x454021d7,0xe3779b90,0xe3779b90,0x81af1549,0x81af1549,0x1fe68f02,0x1fe68f02,0xbe1e08bb,0xbe1e08bb,0x5c558274,0x5c558274,0xfa8cfc2d,0xfa8cfc2d,0x98c475e6,0x98c475e6,0x36fbef9f,0x36fbef9f,0xd5336958,0xd5336958,0x736ae311,0x736ae311,0x11a25cca,0x11a25cca,0xafd9d683,0xafd9d683,0x4e11503c,0x4e11503c,0xec48c9f5,0xec48c9f5,0x8a8043ae,0x8a8043ae,0x28b7bd67,0x28b7bd67,0xc6ef3720});
	//Just precomputed those, might as well

	uint nbRound = k.size();

	for(uint i = 0; i < nbRound; i++){
		if(i%2 == 0)
			addModAddConstValueConstr(model, mk[roundConstants[i]&3], roundConstants[i], k[i], "carryKS"+to_string(i));
		else
			addModAddConstValueConstr(model, mk[(roundConstants[i] >> 11)&3], roundConstants[i], k[i], "carryKS"+to_string(i));
	}
}

GRBModel getXteaModel(uint const nbRound,
					  uint const gamma,
					  GRBEnv & env){

	uint const keySize = 128;

	// Create an empty model
    GRBModel model = GRBModel(env);

    //Key variables
    uint const m = keySize/32;
    vector<vector<GRBVar>> mk(m, vector<GRBVar>(32));
    vector<vector<GRBVar>> k(nbRound, vector<GRBVar>(32));
    for(uint i = 0; i < m; i++){
    	for(uint j = 0; j < 32; j++)
    		mk[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "mk"+to_string(i)+"_"+to_string(j));
    }
    for(uint i = 0; i < nbRound; i++){
    	for(uint j = 0; j < 32; j++)
    		k[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "k"+to_string(i)+"_"+to_string(j));
    }


    //State variables
    vector<vector<GRBVar>> x(nbRound+1, vector<GRBVar>(32));
    vector<vector<GRBVar>> y(nbRound+1, vector<GRBVar>(32));
    for(uint i = 0; i < nbRound+1; i++){
    	for(uint j = 0; j < 32; j++){
    		x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i)+"_"+to_string(j));
    		y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i)+"_"+to_string(j));
    	}
    }
    //Temp vars in the middle of each round
    vector<vector<GRBVar>> yl(nbRound, vector<GRBVar>(32));
    vector<vector<GRBVar>> yr(nbRound, vector<GRBVar>(32));
    vector<vector<GRBVar>> a(nbRound, vector<GRBVar>(32));
    vector<vector<GRBVar>> b(nbRound, vector<GRBVar>(32));
    vector<vector<GRBVar>> c(nbRound, vector<GRBVar>(32));
    for(uint i = 0; i < nbRound; i++){
    	for(uint j = 0; j < 32; j++){
    		yl[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"yl"+to_string(i)+"_"+to_string(j));
    		yr[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"yr"+to_string(i)+"_"+to_string(j));
    		a[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"a"+to_string(i)+"_"+to_string(j));
    		b[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"b"+to_string(i)+"_"+to_string(j));
    		c[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"c"+to_string(i)+"_"+to_string(j));
    	}
    }

    model.update();

    //Key schedule constraints
    addXteaKSConstr(model,mk,k);

    /*
	input x[i] y[i]
	yl[i] = y[i] << 4
	yr[i] = y[i] >> 5
	a[i] = yl[i] ^ yr[i]
	b[i] = a[i] + y[i]
	c[i] = b[i] ^ roundkey[i]
	y[i+1] = x[i] + c[i]
	x[i+1] = y[i]
	output x[i+1] y[i+1]
	*/

    //Round functions
    for(uint r = 0; r < nbRound; r++){

    	//Step1: yl[i] = y[i] << 4; yr[i] = y[i] >> 5; a[i] = yl[i] ^ yr[i]
    	addSHLRXDiffConstr(model,y[r],yl[r],4,gamma);
    	addSHRRXDiffConstr(model,y[r],yr[r],5,gamma);
    	for(uint i = 0; i < 32; i++)
    		addXORConstr(model,yl[r][i],yr[r][i],a[r][i]);

    	//Step2: b[i] = a[i] + y[i]
    	addModAddRXDiffConstr(model,a[r],y[r],b[r],gamma);

    	//Step3: c[i] = b[i] ^ roundkey[i]
    	addRXCstXORConstr(model, b[r], k[r], c[r], gamma);

    	//Step4: y[i+1] = x[i] + c[i]
    	addModAddRXDiffConstr(model,x[r],c[r],y[r+1],gamma);

    	//Step5: x[i+1] = y[i]
    	for(uint i = 0; i < 32; i++)
    		model.addConstr(x[r+1][i] == y[r][i]);
    }

    model.update();

    return model;
}

GRBModel getXteaModelNoKEy(uint const nbRound,
						   uint const gamma,
						   GRBEnv & env){

	// Create an empty model
    GRBModel model = GRBModel(env);

    //State variables
    vector<vector<GRBVar>> x(nbRound+1, vector<GRBVar>(32));
    vector<vector<GRBVar>> y(nbRound+1, vector<GRBVar>(32));
    for(uint i = 0; i < nbRound+1; i++){
    	for(uint j = 0; j < 32; j++){
    		x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i)+"_"+to_string(j));
    		y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i)+"_"+to_string(j));
    	}
    }
    //Temp vars in the middle of each round
    vector<vector<GRBVar>> yl(nbRound, vector<GRBVar>(32));
    vector<vector<GRBVar>> yr(nbRound, vector<GRBVar>(32));
    vector<vector<GRBVar>> a(nbRound, vector<GRBVar>(32));
    vector<vector<GRBVar>> b(nbRound, vector<GRBVar>(32));
    for(uint i = 0; i < nbRound; i++){
    	for(uint j = 0; j < 32; j++){
    		yl[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"yl"+to_string(i)+"_"+to_string(j));
    		yr[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"yr"+to_string(i)+"_"+to_string(j));
    		a[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"a"+to_string(i)+"_"+to_string(j));
    		b[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,"b"+to_string(i)+"_"+to_string(j));
    	}
    }

    model.update();

    /*
	input x[i] y[i]
	yl[i] = y[i] << 4
	yr[i] = y[i] >> 5
	a[i] = yl[i] ^ yr[i]
	b[i] = a[i] + y[i]
	y[i+1] = x[i] + b[i]
	x[i+1] = y[i]
	output x[i+1] y[i+1]
	*/

    //Round functions
    for(uint r = 0; r < nbRound; r++){

    	//Step1: yl[i] = y[i] << 4; yr[i] = y[i] >> 5; a[i] = yl[i] ^ yr[i]
    	addSHLRXDiffConstr(model,y[r],yl[r],4,gamma);
    	addSHRRXDiffConstr(model,y[r],yr[r],5,gamma);
    	for(uint i = 0; i < 32; i++)
    		addXORConstr(model,yl[r][i],yr[r][i],a[r][i]);

    	//Step2: b[i] = a[i] + y[i]
    	addModAddRXDiffConstr(model,a[r],y[r],b[r],gamma);

    	//Step4: y[i+1] = x[i] + b[i]
    	addModAddRXDiffConstr(model,x[r],b[r],y[r+1],gamma);

    	//Step5: x[i+1] = y[i]
    	for(uint i = 0; i < 32; i++)
    		model.addConstr(x[r+1][i] == y[r][i]);
    }

    model.update();

    return model;
}

void searchWeight1RXDiffXtea(uint const nbRound,
							 std::vector<uint> const & cutIndex){
	auto time_start = std::chrono::high_resolution_clock::now();
	
	cout << "---------------------------------------" << endl;
	cout << "XTEA Search for Weight 1 RXDiff with" << endl;
	cout << "nbRound = " << nbRound << endl;
	cout << "cutIndex = {";
	for(auto const tmp : cutIndex)
		cout << tmp << ",";
	cout << "}" << endl;
	cout << "---------------------------------------" << endl;

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);

	uint const keySize = 128;
	uint64_t maxNbKey = (1ULL << cutIndex.size());
	uint m = keySize/32;
	uint n = 32;

	#pragma omp parallel for
	for(uint gamma = 0; gamma < 32; gamma++){

		// cout << "------ For gamma = " << gamma << " ------" << endl;
		uint nbImpossible = 0;
		vector<vector<uint>> listImpossible(64);

		for(uint indexInput = 0; indexInput < 64; indexInput++){
			for(uint indexOutput = 0; indexOutput < 64; indexOutput++){

				auto model = getXteaModel(nbRound, gamma, env);

				//Grab the variables
				vector<GRBVar> plaintext(64);
				vector<GRBVar> ciphertext(64);
				vector<vector<GRBVar>> mk(m, vector<GRBVar>(32));
				for(uint i = 0; i < n; i++){
					plaintext[i] = model.getVarByName("x0_"+to_string(i));
					plaintext[i+n] = model.getVarByName("y0_"+to_string(i));
					ciphertext[i] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i));
					ciphertext[i+n] = model.getVarByName("y"+to_string(nbRound)+"_"+to_string(i));
				}
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < 32; j++)
						mk[i][j] = model.getVarByName("mk"+to_string(i)+"_"+to_string(j));
				}

				//Input diff
				for(uint i = 0; i < indexInput; i++)
					model.addConstr(plaintext[i] == 0);
				model.addConstr(plaintext[indexInput] == 1);
				for(uint i = indexInput+1; i < 64; i++)
					model.addConstr(plaintext[i] == 0);

				//Output diff
				for(uint i = 0; i < indexOutput; i++)
					model.addConstr(ciphertext[i] == 0);
				model.addConstr(ciphertext[indexOutput] == 1);
				for(uint i = indexOutput+1; i < 64; i++)
					model.addConstr(ciphertext[i] == 0);

				//Arbitrary objective to help the solver
				GRBLinExpr obj = 0;
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < n; j++){
						obj += mk[i][j];
					}
				}
				model.setObjective(obj, GRB_MINIMIZE);
				model.set(GRB_IntParam_SolutionLimit, 1);
				//Focus on feasibility
				model.set(GRB_IntParam_MIPFocus,1);

				//Create the callback
				vector<GRBVar> keyvar(n*m);
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < n; j++){
						keyvar[n*i+j] = model.getVarByName("mk"+to_string(i)+"_"+to_string(j));
					}
				}
				CustomCallback cb(keyvar,cutIndex);

				if(cutIndex.size() > 0){
					model.setCallback(&cb);
					model.set(GRB_IntParam_LazyConstraints,1);
				}

				model.update();
				model.optimize();

				if(model.get(GRB_IntAttr_SolCount) > 0){
					// #pragma omp critical
					// {
					// 	cout << "Found key: ";
					// 	for(uint i = 0; i < m; i++){
					// 		uint64_t valk = 0;
					// 		for(uint j = 0; j < 32; j++){
					// 			valk |=  uint(round(mk[i][j].get(GRB_DoubleAttr_X))) << j;
					// 		}
					// 		cout << setfill('0') << setw(8) << hex << valk << " ";
					// 	}
					// 	cout << endl;

					// 	for(uint r = 0; r < nbRound+1; r++){
					// 		if(r < 10) cout << " ";
					// 		cout << "x" << r << ": ";
					// 		uint64_t valxr = 0;
					// 		uint64_t valyr = 0;
					// 		for(uint i = 0; i < 32; i++){
					// 			GRBVar vx = model.getVarByName("x"+to_string(r)+"_"+to_string(i));
					// 			GRBVar vy = model.getVarByName("y"+to_string(r)+"_"+to_string(i));
					// 			valxr |= uint(round(vx.get(GRB_DoubleAttr_X))) << i;
					// 			valyr |= uint(round(vy.get(GRB_DoubleAttr_X))) << i;
					// 		}
					// 		cout << setfill('0') << setw(8) << hex << valxr << " ";
					// 		cout << "y" << r << ": ";
					// 		cout << setfill('0') << setw(8) << hex << valyr << " ";
					// 		cout << endl;
					// 	}
					// }
					continue;
				}
				else if(cutIndex.size() == 0){
					#pragma omp critical
					{	
						// cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						// cout << "gamma = " << gamma << " ";
						// cout << "indexInput = " << indexInput << " indexOutput = " << indexOutput;
						// cout << " --> Impossible differential found" << endl;
						// cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						listImpossible[indexInput].emplace_back(indexOutput);
						nbImpossible++;
					}
				}
				else if(cb.foundKey.size() < maxNbKey){
					#pragma omp critical
					{	
						cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						cout << "indexInput = " << indexInput << " indexOutput = " << indexOutput;
						cout << " --> Impossible differential found" << endl;

						cout << cb.foundKey.size() << " key patterns eliminated" << endl;
						for(auto const & v : cb.foundKey){
							vector<string> tmp(n*m,"*");
							for(auto const & i : cutIndex)
								tmp[i] = to_string(v[i]);
							for(auto const & tt : tmp)
								cout << tt;
							cout << endl;
						}
						cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						listImpossible[indexInput].emplace_back(indexOutput);
						nbImpossible++;
					}
				}
				// else
				// 	cout << "Nope, all keys exhausted" << endl;
			}
		}
		#pragma omp critical
		{
			cout << nbImpossible << " impossible weight (1,1) RX differentials for gamma = " << gamma;
			if(nbImpossible > 0){
				cout << " <<<<<<<<<<<<<<<<<<<<<<<<<";
				cout << endl;

				for(uint indexInput = 0; indexInput < 2*n; indexInput++){
					cout << indexInput << " -> {";
					for(auto const tmp : listImpossible[indexInput])
						cout << tmp << ",";
					cout << "} (" << listImpossible[indexInput].size() << " imp)" << endl;
				}
				cout << "***************************************************" << endl;
			}
			cout << endl;
		}
	}

	auto time_end = std::chrono::high_resolution_clock::now();
	auto total = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	cout << total.count() << " seconds" << endl;

}

void searchWeight1RXDiffXteaNoKey(uint const nbRound){

	auto time_start = std::chrono::high_resolution_clock::now();
	
	cout << "---------------------------------------" << endl;
	cout << "XTEA Search for Weight 1 RXDiff (without key) with" << endl;
	cout << "nbRound = " << nbRound << endl;
	cout << "---------------------------------------" << endl;

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);

	uint n = 32;

	#pragma omp parallel for
	for(uint gamma = 1; gamma < 32; gamma++){

		// cout << "------ For gamma = " << gamma << " ------" << endl;
		uint nbImpossible = 0;

		for(uint indexInput = 0; indexInput < 64; indexInput++){
			for(uint indexOutput = 0; indexOutput < 64; indexOutput++){

				auto model = getXteaModelNoKEy(nbRound, gamma, env);

				//Grab the variables
				vector<GRBVar> plaintext(64);
				vector<GRBVar> ciphertext(64);
				for(uint i = 0; i < n; i++){
					plaintext[i] = model.getVarByName("x0_"+to_string(i));
					plaintext[i+n] = model.getVarByName("y0_"+to_string(i));
					ciphertext[i] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i));
					ciphertext[i+n] = model.getVarByName("y"+to_string(nbRound)+"_"+to_string(i));
				}

				//Input diff
				for(uint i = 0; i < indexInput; i++)
					model.addConstr(plaintext[i] == 0);
				model.addConstr(plaintext[indexInput] == 1);
				for(uint i = indexInput+1; i < 64; i++)
					model.addConstr(plaintext[i] == 0);

				//Output diff
				for(uint i = 0; i < indexOutput; i++)
					model.addConstr(ciphertext[i] == 0);
				model.addConstr(ciphertext[indexOutput] == 1);
				for(uint i = indexOutput+1; i < 64; i++)
					model.addConstr(ciphertext[i] == 0);

				//Arbitrary objective to help the solver
				GRBLinExpr obj = 0;
				for(uint j = 0; j < 32; j++){
					obj += model.getVarByName("x"+to_string(nbRound/2)+"_"+to_string(j));
					obj += model.getVarByName("y"+to_string(nbRound/2)+"_"+to_string(j));
				}
				
				model.setObjective(obj, GRB_MINIMIZE);
				model.set(GRB_IntParam_SolutionLimit, 1);
				//Focus on feasibility
				model.set(GRB_IntParam_MIPFocus,1);

				model.update();
				model.optimize();

				if(model.get(GRB_IntAttr_SolCount) > 0){
					// #pragma omp critical
					// {
					// 	for(uint r = 0; r < nbRound+1; r++){
					// 		if(r < 10) cout << " ";
					// 		cout << "x" << r << ": ";
					// 		uint64_t valxr = 0;
					// 		uint64_t valyr = 0;
					// 		for(uint i = 0; i < 32; i++){
					// 			GRBVar vx = model.getVarByName("x"+to_string(r)+"_"+to_string(i));
					// 			GRBVar vy = model.getVarByName("y"+to_string(r)+"_"+to_string(i));
					// 			valxr |= uint(round(vx.get(GRB_DoubleAttr_X))) << i;
					// 			valyr |= uint(round(vy.get(GRB_DoubleAttr_X))) << i;
					// 		}
					// 		cout << setfill('0') << setw(8) << hex << valxr << " ";
					// 		cout << "y" << r << ": ";
					// 		cout << setfill('0') << setw(8) << hex << valyr << " ";
					// 		cout << endl;
					// 	}
					// }
					continue;
				}
				else{
					#pragma omp critical
					{	
						cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						cout << "gamma = " << gamma << " ";
						cout << "indexInput = " << indexInput << " indexOutput = " << indexOutput;
						cout << " --> Impossible differential found" << endl;
						cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						nbImpossible++;
					}
				}
			}
		}
		#pragma omp critical
		{
			cout << nbImpossible << " impossible weight (1,1) RX differentials for gamma = " << gamma;
			if(nbImpossible > 0)
				cout << "<<<<<<<<<<<<<<<<<<<<<<<<<";
			cout << endl;
		}
	}

	auto time_end = std::chrono::high_resolution_clock::now();
	auto total = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	cout << total.count() << " seconds" << endl;

}
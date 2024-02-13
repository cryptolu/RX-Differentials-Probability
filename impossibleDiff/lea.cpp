#include "lea.hpp"
using namespace std;

void addLeaKSConstr(GRBModel & model,
					std::vector<std::vector<GRBVar>> & mk,
					std::vector<std::vector<std::vector<GRBVar>>> & k){
//Add constraints for the key schedule of LEA
//mk are the variables for the master key (split in m vectors of 32 vars)
//k are the round key variables, split in r*6*32 vars
//These are constraints *in values* not in differences, so constants have an impact

	uint m = mk.size();
	uint nbRound = k.size();

	static const uint delta[] = {0xc3efe9db, 0x44626b02, 0x79e27c8a, 0x78df30ec, 0x715ea49e, 0xc785da0a, 0xe04ef22a, 0xe5c40957};

	//Create variables for the key schedule 
	//We need (r+1)*m*32 T variables
	vector<vector<vector<GRBVar>>> t(nbRound+1, vector<vector<GRBVar>>(m, vector<GRBVar>(32)));
	for(uint r= 0; r < nbRound+1; r++){
		for(uint i = 0; i < m; i++){
			for(uint j = 0; j < 32; j++)
				t[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "t"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
		}
	}

	//Bind the first T vars to the masterkey
	for(uint i = 0; i < m; i++){
		for(uint j = 0; j < 32; j++){
			model.addConstr(t[0][i][j] == mk[i][j]);
		}
	}

	//Key-schedule constraints
	if(m == 4){ //LEA-128
		for(uint i = 0; i < nbRound; i++){

			//State update
			uint cst = CSHL(delta[i%4],i,32);
			vector<GRBVar> out(32);
			//T[0] = T[0]+cst <<< 1
			//To avoid temp variables, this is the same as t[0] >>> 1 = t[0]+cst
			//e.g. (t4,t3,t2,t1,t0) = (x3, x2, x1, x0, x4)
			//<=>  (t0,t4,t3,t2,t1) = (x4, x3, x2, x1, x0)
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][0][(j+1)%32];
			addModAddConstValueConstr(model, t[i][0], cst, out, "carryKS0"+to_string(i));

			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][1][(j+3)%32];
			addModAddConstValueConstr(model, t[i][1], cst, out, "carryKS1"+to_string(i));

			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][2][(j+6)%32];
			addModAddConstValueConstr(model, t[i][2], cst, out, "carryKS2"+to_string(i));

			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][3][(j+11)%32];
			addModAddConstValueConstr(model, t[i][3], cst, out, "carryKS3"+to_string(i));

			//Round key
			for(uint j = 0; j < 32; j++){
				model.addConstr(k[i][0][j] == t[i+1][0][j]);
				model.addConstr(k[i][1][j] == t[i+1][1][j]);
				model.addConstr(k[i][2][j] == t[i+1][2][j]);
				model.addConstr(k[i][3][j] == t[i+1][1][j]);
				model.addConstr(k[i][4][j] == t[i+1][3][j]);
				model.addConstr(k[i][5][j] == t[i+1][1][j]);
			}
		}
	}

	else if(m == 6){
		for(uint i = 0; i < nbRound; i++){

			//State update
			uint cst = CSHL(delta[i%6],i,32);
			vector<GRBVar> out(32);

			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][0][(j+1)%32];
			addModAddConstValueConstr(model, t[i][0], cst, out, "carryKS0"+to_string(i));

			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][1][(j+3)%32];
			addModAddConstValueConstr(model, t[i][1], cst, out, "carryKS1"+to_string(i));

			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][2][(j+6)%32];
			addModAddConstValueConstr(model, t[i][2], cst, out, "carryKS2"+to_string(i));

			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][3][(j+11)%32];
			addModAddConstValueConstr(model, t[i][3], cst, out, "carryKS3"+to_string(i));

			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][4][(j+13)%32];
			addModAddConstValueConstr(model, t[i][4], cst, out, "carryKS4"+to_string(i));

			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][5][(j+17)%32];
			addModAddConstValueConstr(model, t[i][5], cst, out, "carryKS5"+to_string(i));

			//Round key
			for(uint j = 0; j < 32; j++){
				for(uint l = 0; l < 6; l++)
					model.addConstr(k[i][l][j] == t[i+1][l][j]);
			}
		}
	}
	else if(m == 8){
		for(uint i = 0; i < nbRound; i++){

			//State update
			//This case is a bit fucky, as the key schedule only updates 6 words over 8
			//For a much easier time modeling all of this, we keep 8 words for every round, but need to check which ones aren't modified and just add equality constraints

			vector<bool> modifiedWord(8, false); //Ew vector<bool>

			uint cst = CSHL(delta[i%8],i,32);
			vector<GRBVar> out(32);

			uint wordIndex = (6*i)%8;
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][wordIndex][(j+1)%32];
			addModAddConstValueConstr(model, t[i][wordIndex], cst, out, "carryKS0"+to_string(i));
			modifiedWord[wordIndex] = true;

			wordIndex = (6*i+1)%8;
			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][wordIndex][(j+3)%32];
			addModAddConstValueConstr(model, t[i][wordIndex], cst, out, "carryKS1"+to_string(i));
			modifiedWord[wordIndex] = true;

			wordIndex = (6*i+2)%8;
			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][wordIndex][(j+6)%32];
			addModAddConstValueConstr(model, t[i][wordIndex], cst, out, "carryKS2"+to_string(i));
			modifiedWord[wordIndex] = true;

			wordIndex = (6*i+3)%8;
			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][wordIndex][(j+11)%32];
			addModAddConstValueConstr(model, t[i][wordIndex], cst, out, "carryKS3"+to_string(i));
			modifiedWord[wordIndex] = true;

			wordIndex = (6*i+4)%8;
			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][wordIndex][(j+13)%32];
			addModAddConstValueConstr(model, t[i][wordIndex], cst, out, "carryKS4"+to_string(i));
			modifiedWord[wordIndex] = true;

			wordIndex = (6*i+5)%8;
			cst = CSHL(cst,1,32);
			for(uint j = 0; j < 32; j++)
				out[j] = t[i+1][wordIndex][(j+17)%32];
			addModAddConstValueConstr(model, t[i][wordIndex], cst, out, "carryKS5"+to_string(i));
			modifiedWord[wordIndex] = true;

			//Words that are not changed
			for(uint j = 0; j < 8; j++){
				if(!modifiedWord[j]){
					for(uint l = 0; l < 32; l++)
						model.addConstr(t[i+1][j][l] == t[i][j][l]);
				}
			}

			//Round keys
			for(uint j = 0; j < 32; j++){
				for(uint l = 0; l < 6; l++)
					model.addConstr(k[i][l][j] == t[i+1][(6*i+l)%8][j]);
			}

		}
	}
	else{
		cerr << "Error: key size invalid" << endl;
	}
}

GRBModel getLeaModel(uint const keySize,
					 uint const nbRound,
					 uint const gamma,
					 GRBEnv & env){

	/*
	Round function is the following
	x[i+1][0] = CSHL(((x[i][0]^k[i][0]) + (x[i][1]^k[i][1])),9,32)
	x[i+1][1] = CSHR(((x[i][1]^k[i][2]) + (x[i][2]^k[i][3])),5,32)
	x[i+1][2] = CSHR(((x[i][2]^k[i][4]) + (x[i][3]^k[i][5])),3,32)
	x[i+1][3] = x[i][0]
	This is split in the following operations:
	Step1: Constant XOR with the key
	y[i][0] = x[i][0]^k[i][0]
	y[i][1] = x[i][1]^k[i][1]
	y[i][2] = x[i][1]^k[i][2]
	y[i][3] = x[i][2]^k[i][3]
	y[i][4] = x[i][2]^k[i][4]
	y[i][5] = x[i][3]^k[i][5]

	Step2: Modular addition + shift
	x[i+1][0] = CSHL(y[i][0] + y[i][1], 9)
	x[i+1][1] = CSHR(y[i][2] + y[i][3], 5)
	x[i+1][2] = CSHR(y[i][4] + y[i][5], 3)
	x[i+1][3] = x[i][0]
	To avoid creating useless variables, same trick as in the KS, we can use the fact that the relation y = x <<< i is the same as y >>> i = x
	*/

	// Create an empty model
    GRBModel model = GRBModel(env);

    //Key variables
    uint const m = keySize/32;
    vector<vector<GRBVar>> mk(m, vector<GRBVar>(32));
    vector<vector<vector<GRBVar>>> k(nbRound, vector<vector<GRBVar>>(6, vector<GRBVar>(32)));
    for(uint i = 0; i < m; i++){
    	for(uint j = 0; j < 32; j++)
    		mk[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "mk"+to_string(i)+"_"+to_string(j));
    }
    for(uint r = 0; r < nbRound; r++){
    	for(uint i = 0; i < 6; i++){
    		for(uint j = 0; j < 32; j++)
    			k[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "k"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
    	}
    }

    //State variables
    vector<vector<vector<GRBVar>>> x(nbRound+1, vector<vector<GRBVar>>(4, vector<GRBVar>(32)));
    vector<vector<vector<GRBVar>>> y(nbRound, vector<vector<GRBVar>>(6, vector<GRBVar>(32)));
    for(uint r = 0; r < nbRound+1; r++){
    	for(uint i = 0; i < 4; i++){
    		for(uint j = 0; j < 32; j++)
    			x[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
    	}
    }
    for(uint r = 0; r < nbRound; r++){
    	for(uint i = 0; i < 6; i++){
    		for(uint j = 0; j < 32; j++)
    			y[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
    	}
    }

    //Key schedule constraints
    addLeaKSConstr(model,mk,k);

    //Round functions
    for(uint r = 0; r < nbRound; r++){
    	//Step 1
    	addRXCstXORConstr(model, x[r][0], k[r][0], y[r][0], gamma);
    	addRXCstXORConstr(model, x[r][1], k[r][1], y[r][1], gamma);
    	addRXCstXORConstr(model, x[r][1], k[r][2], y[r][2], gamma);
    	addRXCstXORConstr(model, x[r][2], k[r][3], y[r][3], gamma);
    	addRXCstXORConstr(model, x[r][2], k[r][4], y[r][4], gamma);
    	addRXCstXORConstr(model, x[r][3], k[r][5], y[r][5], gamma);

    	//Step 2:
    	vector<GRBVar> out(32); 
    	for(uint i = 0; i < 32; i++)
    		out[i] = x[r+1][0][(i+9)%32]; //out = CSHR(x[i+1][0],9)
    	addModAddRXDiffConstr(model,y[r][0],y[r][1],out,gamma);

    	for(uint i = 0; i < 32; i++)
    		out[i] = x[r+1][1][mod(i-5,32)]; //out = CSHL(x[i+1][1],5)
    	addModAddRXDiffConstr(model,y[r][2],y[r][3],out,gamma);

    	for(uint i = 0; i < 32; i++)
    		out[i] = x[r+1][2][mod(i-3,32)]; //out = CSHL(x[i+1][1],3)
    	addModAddRXDiffConstr(model,y[r][4],y[r][5],out,gamma);

    	for(uint i = 0; i < 32; i++)
    		model.addConstr(x[r+1][3][i] == x[r][0][i]);
    }
    model.update();

    return model;
}

GRBModel getLeaModelNoKey(uint const nbRound,
						  uint const gamma,
						  GRBEnv & env){

	/*
	Round function is the following
	x[i+1][0] = CSHL(((x[i][0]^k[i][0]) + (x[i][1]^k[i][1])),9,32)
	x[i+1][1] = CSHR(((x[i][1]^k[i][2]) + (x[i][2]^k[i][3])),5,32)
	x[i+1][2] = CSHR(((x[i][2]^k[i][4]) + (x[i][3]^k[i][5])),3,32)
	x[i+1][3] = x[i][0]

	here we ignore the key, so it's just

	Step2: Modular addition + shift
	x[i+1][0] = CSHL(x[i][0] + x[i][1], 9)
	x[i+1][1] = CSHR(x[i][1] + x[i][2], 5)
	x[i+1][2] = CSHR(x[i][2] + x[i][3], 3)
	x[i+1][3] = x[i][0]
	To avoid creating useless variables, same trick as in the KS, we can use the fact that the relation y = x <<< i is the same as y >>> i = x
	*/

	// Create an empty model
    GRBModel model = GRBModel(env);

    //State variables
    vector<vector<vector<GRBVar>>> x(nbRound+1, vector<vector<GRBVar>>(4, vector<GRBVar>(32)));
    for(uint r = 0; r < nbRound+1; r++){
    	for(uint i = 0; i < 4; i++){
    		for(uint j = 0; j < 32; j++)
    			x[r][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
    	}
    }

    //Round functions
    for(uint r = 0; r < nbRound; r++){
    	//Step 2:
    	vector<GRBVar> out(32); 
    	for(uint i = 0; i < 32; i++)
    		out[i] = x[r+1][0][(i+9)%32]; //out = CSHR(x[i+1][0],9)
    	addModAddRXDiffConstr(model,x[r][0],x[r][1],out,gamma);

    	for(uint i = 0; i < 32; i++)
    		out[i] = x[r+1][1][mod(i-5,32)]; //out = CSHL(x[i+1][1],5)
    	addModAddRXDiffConstr(model,x[r][1],x[r][2],out,gamma);

    	for(uint i = 0; i < 32; i++)
    		out[i] = x[r+1][2][mod(i-3,32)]; //out = CSHL(x[i+1][1],3)
    	addModAddRXDiffConstr(model,x[r][2],x[r][3],out,gamma);

    	for(uint i = 0; i < 32; i++)
    		model.addConstr(x[r+1][3][i] == x[r][0][i]);
    }
    model.update();

    return model;
}

void searchWeight1RXDiffLea(uint const keySize,
							uint const nbRound,
							std::vector<uint> const & cutIndex){
	auto time_start = std::chrono::high_resolution_clock::now();
	
	cout << "---------------------------------------" << endl;
	cout << "LEA Search for Weight 1 RXDiff with" << endl;
	cout << "keySize = " << keySize << " nbRound = " << nbRound << endl;
	cout << "cutIndex = {";
	for(auto const tmp : cutIndex)
		cout << tmp << ",";
	cout << "}" << endl;
	cout << "---------------------------------------" << endl;

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);

	uint64_t maxNbKey = (1ULL << cutIndex.size());
	uint m = keySize/32;
	uint n = 32;

	#pragma omp parallel for
	for(uint gamma = 1; gamma < 32; gamma++){

		// cout << "------ For gamma = " << gamma << " ------" << endl;
		uint nbImpossible = 0;

		for(uint indexInput = 0; indexInput < 128; indexInput++){
			for(uint indexOutput = 0; indexOutput < 128; indexOutput++){

				auto model = getLeaModel(keySize, nbRound, gamma, env);

				//Grab the variables
				vector<GRBVar> plaintext(128);
				vector<GRBVar> ciphertext(128);
				vector<vector<GRBVar>> mk(m, vector<GRBVar>(32));
				uint ctr = 0;
				for(uint i = 0; i < 4; i++){
					for(uint j = 0; j < 32; j++){
						plaintext[ctr] = model.getVarByName("x0_"+to_string(i)+"_"+to_string(j));
						ciphertext[ctr] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i)+"_"+to_string(j));
						ctr++;
					}
				}
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < 32; j++)
						mk[i][j] = model.getVarByName("mk"+to_string(i)+"_"+to_string(j));
				}

				//Input diff
				for(uint i = 0; i < indexInput; i++)
					model.addConstr(plaintext[i] == 0);
				model.addConstr(plaintext[indexInput] == 1);
				for(uint i = indexInput+1; i < 128; i++)
					model.addConstr(plaintext[i] == 0);

				//Output diff
				for(uint i = 0; i < indexOutput; i++)
					model.addConstr(ciphertext[i] == 0);
				model.addConstr(ciphertext[indexOutput] == 1);
				for(uint i = indexOutput+1; i < 128; i++)
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
					// cout << "Found key: ";
					// for(uint i = 0; i < m; i++){
					// 	uint64_t valk = 0;
					// 	for(uint j = 0; j < 32; j++){
					// 		valk |=  uint(round(mk[i][j].get(GRB_DoubleAttr_X))) << j;
					// 	}
					// 	cout << setfill('0') << setw(8) << hex << valk << " ";
					// }
					// cout << endl;

					// for(uint r = 0; r < nbRound+1; r++){
					// 	if(r < 10) cout << " ";
					// 	cout << "x" << r << ": ";
					// 	for(uint i = 0; i < 4; i++){
					// 		uint64_t valxi = 0;
					// 		for(uint j = 0; j < 32; j++){
					// 			GRBVar v = model.getVarByName("x"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
					// 			valxi |= uint(round(v.get(GRB_DoubleAttr_X))) << j;
					// 		}
					// 		cout << setfill('0') << setw(8) << hex << valxi << " ";
					// 	}
					// 	cout << endl;
					// }
					continue;
				}
				else if(cutIndex.size() == 0){
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
			if(nbImpossible > 0)
				cout << "<<<<<<<<<<<<<<<<<<<<<<<<<";
			cout << endl;
		}
	}

	auto time_end = std::chrono::high_resolution_clock::now();
	auto total = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	cout << total.count() << " seconds" << endl;
}

void searchWeight1RXDiffLeaNoKey(uint const nbRound){
	auto time_start = std::chrono::high_resolution_clock::now();
	
	cout << "---------------------------------------" << endl;
	cout << "LEA Search for Weight 1 RXDiff (without key) with" << endl;
	cout << "nbRound = " << nbRound << endl;
	cout << "---------------------------------------" << endl;

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);

	#pragma omp parallel for
	for(uint gamma = 1; gamma < 32; gamma++){

		// cout << "------ For gamma = " << gamma << " ------" << endl;
		uint nbImpossible = 0;

		for(uint indexInput = 0; indexInput < 128; indexInput++){
			for(uint indexOutput = 0; indexOutput < 128; indexOutput++){

				auto model = getLeaModelNoKey(nbRound, gamma, env);

				//Grab the variables
				vector<GRBVar> plaintext(128);
				vector<GRBVar> ciphertext(128);
				uint ctr = 0;
				for(uint i = 0; i < 4; i++){
					for(uint j = 0; j < 32; j++){
						plaintext[ctr] = model.getVarByName("x0_"+to_string(i)+"_"+to_string(j));
						ciphertext[ctr] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i)+"_"+to_string(j));
						ctr++;
					}
				}

				//Input diff
				for(uint i = 0; i < indexInput; i++)
					model.addConstr(plaintext[i] == 0);
				model.addConstr(plaintext[indexInput] == 1);
				for(uint i = indexInput+1; i < 128; i++)
					model.addConstr(plaintext[i] == 0);

				//Output diff
				for(uint i = 0; i < indexOutput; i++)
					model.addConstr(ciphertext[i] == 0);
				model.addConstr(ciphertext[indexOutput] == 1);
				for(uint i = indexOutput+1; i < 128; i++)
					model.addConstr(ciphertext[i] == 0);

				//Arbitrary objective to help the solver
				GRBLinExpr obj = 0;
				for(uint i = 0; i < 4; i++){
					for(uint j = 0; j < 32; j++){
						obj += model.getVarByName("x"+to_string(nbRound/2)+"_"+to_string(i)+"_"+to_string(j));
						ctr++;
					}
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
					// 		for(uint i = 0; i < 4; i++){
					// 			uint64_t valxi = 0;
					// 			for(uint j = 0; j < 32; j++){
					// 				GRBVar v = model.getVarByName("x"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
					// 				valxi |= uint(round(v.get(GRB_DoubleAttr_X))) << j;
					// 			}
					// 			cout << setfill('0') << setw(8) << hex << valxi << " ";
					// 		}
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
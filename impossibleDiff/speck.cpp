#include "speck.hpp"
using namespace std;

void addSpeckKSConstr(GRBModel & model,
					  std::vector<std::vector<GRBVar>> & mk,
					  std::vector<std::vector<GRBVar>> & k,
					  uint const alpha,
					  uint const beta){
//Add constraints for the key schedule of Speck
//mk are the variables for the master key (split in m vectors)
//k are the variables for the round keys (split in r vectors)
//n (word size) is obtained from mk[0].size()
//These are constraints *in values* not in differences, so constants have an impact

	uint m = mk.size();
	uint n = mk[0].size();
	uint r = k.size();

	//Create variables for the key schedule 
	//We need the l[i] variables from the specification, but also some additional temp variables
	vector<vector<GRBVar>> l(r+m-2, vector<GRBVar>(n)); //l variables
	vector<vector<GRBVar>> tmpKS(r-1, vector<GRBVar>(n)); //tmp variables for the computation of l[i]

	for(uint i = 0; i < r+m-2; i++){
		for(uint j = 0; j < n; j++)
			l[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "l"+to_string(i)+"_"+to_string(j));
	}
	for(uint i = 0; i < r-1; i++){
		for(uint j = 0; j < n; j++)
			tmpKS[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "tmpKS"+to_string(i)+"_"+to_string(j));
	}
	model.update();

	//First, constraints to bind k[0] and the first few l[i] to the master key
	for(uint j = 0; j < n; j++)
		model.addConstr(k[0][j] == mk[0][j]);
	for(uint i = 0; i < m-1; i++){
		for(uint j = 0; j < n; j++){
			model.addConstr(l[i][j] == mk[i+1][j]);
		}
	}

	//Key schedule constraints
	for(uint64_t i = 0; i < r-1; i++){
		vector<GRBVar> Sali(n); //l[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Sali[j] = l[i][mod(j+alpha,n)];

		vector<GRBVar> Sbki(n); //k[i] rotated to the left by beta
		for(uint j = 0; j < n; j++)
			Sbki[j] = k[i][mod(j-beta,n)];

		//l[i+m-1] = (k[i] + Sali) ^ i
		//First, tmp = k[i] + Sali
		addModAddValueConstr(model, k[i], Sali, tmpKS[i], "carryKS"+to_string(i));

		//l[i+m-1] = tmp ^ i 
		for(uint j = 0; j < n; j++){
			if(((i >> j)&1) == 0)
				model.addConstr(l[i+m-1][j] == tmpKS[i][j]);
			else
				model.addConstr(l[i+m-1][j] == 1 - tmpKS[i][j]);
		}

		//k[i+1] = Sbki ^ l[i+m-1]
		for(uint j = 0; j < n; j++)
			addXORConstr(model, Sbki[j], l[i+m-1][j], k[i+1][j]);
	}
}

void addSpeckKSRelatedKeyConstr(GRBModel & model,
								std::vector<std::vector<GRBVar>> & mk,
								std::vector<std::vector<GRBVar>> & k,
								uint const alpha,
								uint const beta,
								uint const gamma){
//Add constraints for the key schedule of Speck
//mk are the variables for the master key (split in m vectors)
//k are the variables for the round keys (split in r vectors)
//n (word size) is obtained from mk[0].size()
//These are constraints *in values* not in differences, so constants have an impact

	uint m = mk.size();
	uint n = mk[0].size();
	uint r = k.size();

	//Create variables for the key schedule 
	//We need the l[i] variables from the specification, but also some additional temp variables
	vector<vector<GRBVar>> l(r+m-2, vector<GRBVar>(n)); //l variables
	vector<vector<GRBVar>> tmpKS(r-1, vector<GRBVar>(n)); //tmp variables for the computation of l[i]

	for(uint i = 0; i < r+m-2; i++){
		for(uint j = 0; j < n; j++)
			l[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "l"+to_string(i)+"_"+to_string(j));
	}
	for(uint i = 0; i < r-1; i++){
		for(uint j = 0; j < n; j++)
			tmpKS[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "tmpKS"+to_string(i)+"_"+to_string(j));
	}
	model.update();

	//First, constraints to bind k[0] and the first few l[i] to the master key
	for(uint j = 0; j < n; j++)
		model.addConstr(k[0][j] == mk[0][j]);
	for(uint i = 0; i < m-1; i++){
		for(uint j = 0; j < n; j++){
			model.addConstr(l[i][j] == mk[i+1][j]);
		}
	}

	//Key schedule constraints
	for(uint64_t i = 0; i < r-1; i++){
		vector<GRBVar> Sali(n); //l[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Sali[j] = l[i][mod(j+alpha,n)];

		vector<GRBVar> Sbki(n); //k[i] rotated to the left by beta
		for(uint j = 0; j < n; j++)
			Sbki[j] = k[i][mod(j-beta,n)];

		//Write the round constant (i) in binary
		vector<uint> roundcst(n,0);
		for(uint j = 0; j < n; j++){
			if(((i >> j)&1) == 1)
				roundcst[j] = 1;
		}

		//l[i+m-1] = (k[i] + Sali) ^ i
		//First, tmp = k[i] + Sali
		addModAddRXDiffConstr(model, k[i], Sali, tmpKS[i], gamma);

		//l[i+m-1] = tmp ^ i 
		addRXCstXORConstr(model, tmpKS[i], roundcst, l[i+m-1], gamma);

		//k[i+1] = Sbki ^ l[i+m-1]
		for(uint j = 0; j < n; j++)
			addXORConstr(model, Sbki[j], l[i+m-1][j], k[i+1][j]);
	}
}


GRBModel
getSpeckModel(uint const n,
			  uint const m, 
			  uint const nbRound,
			  uint const alpha, 
			  uint const beta, 
			  uint const gamma,
			  GRBEnv & env){
//Create a model for single-key RX-diff for Speck 2n/mn over nbRounds rounds, using alpha/beta for the round function. env is only here for better console login from gurobi
//Returns a pair with the model, and a vector of length 3 containing the x vars, y vars, and mk vars, in this order

	// Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    //Key variables
    vector<vector<GRBVar>> mk(m, vector<GRBVar>(n)); //Special variables for the master key, for easier use
    vector<vector<GRBVar>> k(nbRound, vector<GRBVar>(n));
	for(uint i = 0; i < nbRound; i++){
		for(uint j = 0; j < n; j++)
			k[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "k"+to_string(i)+"_"+to_string(j));
	}
	for(uint i = 0; i < m; i++){
		for(uint j = 0; j < n; j++)
			mk[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "mk"+to_string(i)+"_"+to_string(j));
	}

	/*
	xi      yi
	|       |
	S-alpha |
	|       |
	+-------|
	zi      |
	|       Sbeta
	^-ki    |
	|-------^
	|       |
	xi+1    yi+1
	*/

	//State variables
	vector<vector<GRBVar>> x(nbRound+1, vector<GRBVar>(n));
	vector<vector<GRBVar>> y(nbRound+1, vector<GRBVar>(n));
	vector<vector<GRBVar>> z(nbRound, vector<GRBVar>(n));
	for(uint i = 0; i < nbRound+1; i++){
		for(uint j = 0; j < n; j++){
			x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i)+"_"+to_string(j));
			y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i)+"_"+to_string(j));
		}
	}
	for(uint i = 0; i < nbRound; i++){
		for(uint j = 0; j < n; j++){
			z[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z"+to_string(i)+"_"+to_string(j));
		}
	}


	//Key schedule constraints
	addSpeckKSConstr(model,mk,k,alpha,beta);

	//Round functions
	for(uint i = 0; i < nbRound; i++){
		vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Saxi[j] = x[i][mod(j+alpha,n)];

		vector<GRBVar> Sbyi(n); //k[i] rotated to the left by beta
		for(uint j = 0; j < n; j++)
			Sbyi[j] = y[i][mod(j-beta,n)];

		//zi = Saxi + yi
		addModAddRXDiffConstr(model,Saxi,y[i],z[i],gamma);

		//x[i+1] = z[i] ^ k[i]
		addRXCstXORConstr(model, z[i], k[i], x[i+1], gamma);

		//y[i+1] = Sbyi ^ x[i+1]
		for(uint j = 0; j < n; j++)
			addXORConstr(model, Sbyi[j], x[i+1][j], y[i+1][j]);
	}

	return model;
}

GRBModel
getSpeckModelNoKey(uint const n,
				   uint const nbRound,
				   uint const alpha, 
				   uint const beta, 
				   uint const gamma,
				   GRBEnv & env){
//Create a model for single-key RX-diff for Speck 2n/mn over nbRounds rounds, using alpha/beta for the round function. env is only here for better console login from gurobi

	// Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables

	/*
	xi      yi
	|       |
	S-alpha |
	|       |
	+-------|
	|       Sbeta
	|       |
	|-------^
	|       |
	xi+1    yi+1
	*/

	//State variables
	vector<vector<GRBVar>> x(nbRound+1, vector<GRBVar>(n));
	vector<vector<GRBVar>> y(nbRound+1, vector<GRBVar>(n));
	for(uint i = 0; i < nbRound+1; i++){
		for(uint j = 0; j < n; j++){
			x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i)+"_"+to_string(j));
			y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i)+"_"+to_string(j));
		}
	}

	//Round functions
	for(uint i = 0; i < nbRound; i++){
		vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Saxi[j] = x[i][mod(j+alpha,n)];

		vector<GRBVar> Sbyi(n); //k[i] rotated to the left by beta
		for(uint j = 0; j < n; j++)
			Sbyi[j] = y[i][mod(j-beta,n)];

		//x[i+1] = Saxi + yi
		addModAddRXDiffConstr(model,Saxi,y[i],x[i+1],gamma);

		//y[i+1] = Sbyi ^ x[i+1]
		for(uint j = 0; j < n; j++)
			addXORConstr(model, Sbyi[j], x[i+1][j], y[i+1][j]);
	}
	model.update();

	return model;
}

GRBModel
getSpeckRelatedKeyModel(uint const n,
						uint const m, 
						uint const nbRound,
						uint const alpha, 
						uint const beta, 
						uint const gamma,
						GRBEnv & env){
//Create a model for single-key RX-diff for Speck 2n/mn over nbRounds rounds, using alpha/beta for the round function. env is only here for better console login from gurobi
//Returns a pair with the model, and a vector of length 3 containing the x vars, y vars, and mk vars, in this order

	// Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    //Key variables
    vector<vector<GRBVar>> mk(m, vector<GRBVar>(n)); //Special variables for the master key, for easier use
    vector<vector<GRBVar>> k(nbRound, vector<GRBVar>(n));
	for(uint i = 0; i < nbRound; i++){
		for(uint j = 0; j < n; j++)
			k[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "k"+to_string(i)+"_"+to_string(j));
	}
	for(uint i = 0; i < m; i++){
		for(uint j = 0; j < n; j++)
			mk[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "mk"+to_string(i)+"_"+to_string(j));
	}

	/*
	xi      yi
	|       |
	S-alpha |
	|       |
	+-------|
	zi      |
	|       Sbeta
	^-ki    |
	|-------^
	|       |
	xi+1    yi+1
	*/

	//State variables
	vector<vector<GRBVar>> x(nbRound+1, vector<GRBVar>(n));
	vector<vector<GRBVar>> y(nbRound+1, vector<GRBVar>(n));
	vector<vector<GRBVar>> z(nbRound, vector<GRBVar>(n));
	for(uint i = 0; i < nbRound+1; i++){
		for(uint j = 0; j < n; j++){
			x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i)+"_"+to_string(j));
			y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i)+"_"+to_string(j));
		}
	}
	for(uint i = 0; i < nbRound; i++){
		for(uint j = 0; j < n; j++){
			z[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z"+to_string(i)+"_"+to_string(j));
		}
	}


	//Key schedule constraints
	addSpeckKSRelatedKeyConstr(model,mk,k,alpha,beta,gamma);

	//Round functions
	for(uint i = 0; i < nbRound; i++){
		vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Saxi[j] = x[i][mod(j+alpha,n)];

		vector<GRBVar> Sbyi(n); //k[i] rotated to the left by beta
		for(uint j = 0; j < n; j++)
			Sbyi[j] = y[i][mod(j-beta,n)];

		//zi = Saxi + yi
		addModAddRXDiffConstr(model,Saxi,y[i],z[i],gamma);

		//x[i+1] = z[i] ^ k[i]
		//y[i+1] = Sbyi ^ x[i+1]
		for(uint j = 0; j < n; j++){
			addXORConstr(model, z[i][j], k[i][j], x[i+1][j]);
			addXORConstr(model, Sbyi[j], x[i+1][j], y[i+1][j]);
		}
	}

	return model;
}


void searchWeight1RXDiffSpeck(uint const n,
							  uint const m,
							  uint const nbRound,
							  std::vector<uint> const & cutIndex){

	auto time_start = std::chrono::high_resolution_clock::now();


	cout << "---------------------------------------" << endl;
	cout << "Speck Search for Weight 1 RXDiff with" << endl;
	cout << "n = " << n << " m = " << m << " nbRound = " << nbRound << endl;
	cout << "cutIndex = {";
	for(auto const tmp : cutIndex)
		cout << tmp << ",";
	cout << "}" << endl;
	cout << "---------------------------------------" << endl;

	//Modelize and solve
	uint alpha = 8;
	uint beta = 3;
	if(n == 16){
		alpha = 7;
		beta = 2;
	}
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);

	uint64_t maxNbKey = (1ULL << cutIndex.size());

	#pragma omp parallel for
	for(uint gamma = 1; gamma < n; gamma++){
		// cout << "------ For gamma = " << gamma << " ------" << endl;
		uint nbImpossible = 0;

		for(uint indexInput = 0; indexInput < 2*n; indexInput++){
			for(uint indexOutput = 0; indexOutput < 2*n; indexOutput++){

				// cout << "indexInput = " << indexInput << " indexOutput = " << indexOutput << endl;

				auto model = getSpeckModel(n,m,nbRound,alpha,beta,gamma,env);

				//Grab the variables
				vector<GRBVar> plaintext(2*n);
				vector<GRBVar> ciphertext(2*n);
				for(uint i = 0; i < n; i++){
					plaintext[i] = model.getVarByName("x0_"+to_string(i));
					plaintext[i+n] = model.getVarByName("y0_"+to_string(i));
					ciphertext[i] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i));
					ciphertext[i+n] = model.getVarByName("y"+to_string(nbRound)+"_"+to_string(i));
				}

				vector<vector<GRBVar>> mk(m, vector<GRBVar>(n));
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < n; j++){
						mk[i][j] = model.getVarByName("mk"+to_string(i)+"_"+to_string(j));
					}
				}

				//Input diff
				for(uint i = 0; i < indexInput; i++)
					model.addConstr(plaintext[i] == 0);
				model.addConstr(plaintext[indexInput] == 1);
				for(uint i = indexInput+1; i < 2*n; i++)
					model.addConstr(plaintext[i] == 0);

				//Output diff
				for(uint i = 0; i < indexOutput; i++)
					model.addConstr(ciphertext[i] == 0);
				model.addConstr(ciphertext[indexOutput] == 1);
				for(uint i = indexOutput+1; i < 2*n; i++)
					model.addConstr(ciphertext[i] == 0);

				//Arbitrary objective to help the solver
				GRBLinExpr obj = 0;
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < n; j++){
						obj += mk[i][j];
					}
				}

				//Create the callback
				vector<GRBVar> keyvar(n*m);
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < n; j++){
						keyvar[n*i+j] = model.getVarByName("mk"+to_string(i)+"_"+to_string(j));
					}
				}
				model.update();
				CustomCallback cb(keyvar,cutIndex);
				
				if(cutIndex.size() > 0){
					model.setCallback(&cb);
					model.set(GRB_IntParam_LazyConstraints,1);
				}

				model.setObjective(obj, GRB_MINIMIZE);
				model.set(GRB_IntParam_SolutionLimit, 1);
				//Focus on feasibility
				model.set(GRB_IntParam_MIPFocus,1);

				model.update();
				model.optimize();

				if(model.get(GRB_IntAttr_SolCount) > 0){
					//Keep this for checking with prints if necessary
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
						cout << "gamma = " << gamma << " ";
						cout << "indexInput = " << indexInput << " indexOutput = " << indexOutput << endl;
						cout << "Impossible differential found" << endl;

						cout << dec;

						cout << cb.foundKey.size() << " key patterns eliminated" << endl;
						for(auto const & v : cb.foundKey){
							vector<string> tmp(n*m,"*");
							for(auto const & i : cutIndex)
								tmp[i] = to_string(v[i]);
							for(auto const & tt : tmp)
								cout << tt;
							cout << endl;
						}
						nbImpossible++;
						cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
					}
				}
				// else{
				// 	cout << "Nope, all keys exhausted" << endl;
				// }
			}
		}

		#pragma omp critical
		{
			// cout << nbPossible << " possible weight (1,1) RX differentials" << endl;
			cout << nbImpossible << " impossible weight (1,1) RX differentials";
			if(nbImpossible > 0)
				cout << "<<<<<<<<<<<<<<<<<<<<<<<<<";
			cout << endl;
		}
	}

	auto time_end = std::chrono::high_resolution_clock::now();
	auto total = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	cout << total.count() << " seconds" << endl;
}

void searchWeight1RXDiffSpeckNoKey(uint const n,
								   uint const nbRound){

	auto time_start = std::chrono::high_resolution_clock::now();


	cout << "---------------------------------------" << endl;
	cout << "Speck Search for Weight 1 RXDiff (without key) with" << endl;
	cout << "n = " << n << " nbRound = " << nbRound << endl;
	cout << "---------------------------------------" << endl;

	//Modelize and solve
	uint alpha = 8;
	uint beta = 3;
	if(n == 16){
		alpha = 7;
		beta = 2;
	}
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);

	#pragma omp parallel for
	for(uint gamma = 0; gamma < n; gamma++){
		// cout << "------ For gamma = " << gamma << " ------" << endl;
		uint nbImpossible = 0;
		vector<vector<uint>> listImpossible(2*n);

		for(uint indexInput = 0; indexInput < 2*n; indexInput++){
			for(uint indexOutput = 0; indexOutput < 2*n; indexOutput++){

				// cout << "indexInput = " << indexInput << " indexOutput = " << indexOutput << endl;

				auto model = getSpeckModelNoKey(n,nbRound,alpha,beta,gamma,env);

				//Grab the variables
				vector<GRBVar> plaintext(2*n);
				vector<GRBVar> ciphertext(2*n);
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
				for(uint i = indexInput+1; i < 2*n; i++)
					model.addConstr(plaintext[i] == 0);

				//Output diff
				for(uint i = 0; i < indexOutput; i++)
					model.addConstr(ciphertext[i] == 0);
				model.addConstr(ciphertext[indexOutput] == 1);
				for(uint i = indexOutput+1; i < 2*n; i++)
					model.addConstr(ciphertext[i] == 0);

				//Arbitrary objective to help the solver
				GRBLinExpr obj = 0;
				for(uint j = 0; j < n; j++){
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
					//Keep this for checking with prints if necessary
					continue;
				}
				else{
					#pragma omp critical
					{	
						// cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						// cout << "gamma = " << gamma << " ";
						// cout << "indexInput = " << indexInput << " indexOutput = " << indexOutput;
						// cout << " --> Impossible differential found" << endl;
						// cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						nbImpossible++;
						listImpossible[indexInput].emplace_back(indexOutput);
					}
				}
			}
		}

		#pragma omp critical
		{
			// cout << nbPossible << " possible weight (1,1) RX differentials" << endl;
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


bool existTrailSpeckNoKey(uint const n,
						  uint const nbRound,
						  uint const gamma,
						  std::vector<uint> const & inputDiff,
						  std::vector<uint> const & outputDiff,
						  GRBEnv & env){

	uint alpha = 8;
	uint beta = 3;
	if(n == 16){
		alpha = 7;
		beta = 2;
	}

	auto model = getSpeckModelNoKey(n,nbRound,alpha,beta,gamma,env);

	//Grab the variables
	vector<GRBVar> plaintext(2*n);
	vector<GRBVar> ciphertext(2*n);
	for(uint i = 0; i < n; i++){
		plaintext[i] = model.getVarByName("x0_"+to_string(i));
		plaintext[i+n] = model.getVarByName("y0_"+to_string(i));
		ciphertext[i] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i));
		ciphertext[i+n] = model.getVarByName("y"+to_string(nbRound)+"_"+to_string(i));
	}

	for(uint i = 0; i < 2*n; i++){
		model.addConstr(plaintext[i] == inputDiff[i]);
		model.addConstr(ciphertext[i] == outputDiff[i]);
	}
	//Arbitrary objective to help the solver
	GRBLinExpr obj = 0;
	for(uint i = 0; i < n; i++){
		obj += model.getVarByName("x"+to_string(nbRound/2)+"_"+to_string(i));
		obj += model.getVarByName("y"+to_string(nbRound/2)+"_"+to_string(i));
	}

	model.setObjective(obj, GRB_MINIMIZE);
	model.set(GRB_IntParam_SolutionLimit, 1);
	//Focus on feasibility
	model.set(GRB_IntParam_MIPFocus,1);

	model.update();
	model.optimize();
	return model.get(GRB_IntAttr_SolCount) > 0;

}

void searchWeight1RXDiffSpeckRelatedKey(uint const n,
										uint const m,
										uint const nbRound){

	auto time_start = std::chrono::high_resolution_clock::now();


	cout << "---------------------------------------" << endl;
	cout << "Speck Search for Weight 1 RXDiff Related Key with" << endl;
	cout << "n = " << n << " m = " << m << " nbRound = " << nbRound << endl;
	cout << "---------------------------------------" << endl;

	//Modelize and solve
	uint alpha = 8;
	uint beta = 3;
	if(n == 16){
		alpha = 7;
		beta = 2;
	}
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);

	#pragma omp parallel for
	for(uint gamma = 1; gamma < n; gamma++){
		for(uint indexKey = 0; indexKey < n*m; indexKey++){
		// cout << "------ For gamma = " << gamma << " ------" << endl;
		uint nbImpossible = 0;
		vector<vector<uint>> listImpossible(2*n);

		for(uint indexInput = 0; indexInput < 2*n; indexInput++){
			for(uint indexOutput = 0; indexOutput < 2*n; indexOutput++){

				// cout << "indexInput = " << indexInput << " indexOutput = " << indexOutput << endl;

				auto model = getSpeckRelatedKeyModel(n,m,nbRound,alpha,beta,gamma,env);

				//Grab the variables
				vector<GRBVar> plaintext(2*n);
				vector<GRBVar> ciphertext(2*n);
				for(uint i = 0; i < n; i++){
					plaintext[i] = model.getVarByName("x0_"+to_string(i));
					plaintext[i+n] = model.getVarByName("y0_"+to_string(i));
					ciphertext[i] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i));
					ciphertext[i+n] = model.getVarByName("y"+to_string(nbRound)+"_"+to_string(i));
				}

				vector<vector<GRBVar>> mk(m, vector<GRBVar>(n));
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < n; j++){
						mk[i][j] = model.getVarByName("mk"+to_string(i)+"_"+to_string(j));
					}
				}

				//Input diff
				for(uint i = 0; i < indexInput; i++)
					model.addConstr(plaintext[i] == 0);
				model.addConstr(plaintext[indexInput] == 1);
				for(uint i = indexInput+1; i < 2*n; i++)
					model.addConstr(plaintext[i] == 0);

				//Output diff
				for(uint i = 0; i < indexOutput; i++)
					model.addConstr(ciphertext[i] == 0);
				model.addConstr(ciphertext[indexOutput] == 1);
				for(uint i = indexOutput+1; i < 2*n; i++)
					model.addConstr(ciphertext[i] == 0);

				//Arbitrary objective to help the solver
				GRBLinExpr obj = 0;
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < n; j++){
						obj += mk[i][j];
					}
				}

				//Key diff
				vector<GRBVar> keyvar(n*m);
				for(uint i = 0; i < m; i++){
					for(uint j = 0; j < n; j++){
						keyvar[n*i+j] = model.getVarByName("mk"+to_string(i)+"_"+to_string(j));
					}
				}
				for(uint i = 0; i < indexKey; i++)
					model.addConstr(keyvar[i] == 0);
				model.addConstr(keyvar[indexKey] == 1);
				for(uint i = indexKey+1; i < n*m; i++)
					model.addConstr(keyvar[i] == 0);
				model.update();
				

				model.setObjective(obj, GRB_MINIMIZE);
				model.set(GRB_IntParam_SolutionLimit, 1);
				//Focus on feasibility
				model.set(GRB_IntParam_MIPFocus,1);

				model.update();
				model.optimize();

				if(model.get(GRB_IntAttr_SolCount) > 0){
					//Keep this for checking with prints if necessary
					continue;
				}
				else{
					#pragma omp critical
					{	
						// cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						// cout << "gamma = " << gamma << " ";
						// cout << "indexInput = " << indexInput << " indexOutput = " << indexOutput;
						// cout << " --> Impossible differential found" << endl;
						// cout << "!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
						nbImpossible++;
						listImpossible[indexInput].emplace_back(indexOutput);
					}
				}
			}
		}

		#pragma omp critical
		{
			cout << nbImpossible << " impossible weight (1,1) RX differentials for gamma = " << gamma << " and indexKey = " << indexKey;
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
	}}

	auto time_end = std::chrono::high_resolution_clock::now();
	auto total = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	cout << total.count() << " seconds" << endl;
}

bool existTrailSpeckRelatedKey(uint const n,
							   uint const m,
							   uint const nbRound,
							   uint const gamma,
							   std::vector<uint> const & inputDiff,
							   std::vector<uint> const & outputDiff,
							   std::vector<uint> const & keyDiff,
							   GRBEnv & env){

	//Modelize and solve
	uint alpha = 8;
	uint beta = 3;
	if(n == 16){
		alpha = 7;
		beta = 2;
	}

	auto model = getSpeckRelatedKeyModel(n,m,nbRound,alpha,beta,gamma,env);

	//Grab the variables
	vector<GRBVar> plaintext(2*n);
	vector<GRBVar> ciphertext(2*n);
	for(uint i = 0; i < n; i++){
		plaintext[i] = model.getVarByName("x0_"+to_string(i));
		plaintext[i+n] = model.getVarByName("y0_"+to_string(i));
		ciphertext[i] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i));
		ciphertext[i+n] = model.getVarByName("y"+to_string(nbRound)+"_"+to_string(i));
	}
	vector<GRBVar> keyvar(n*m);
	for(uint i = 0; i < m; i++){
		for(uint j = 0; j < n; j++){
			keyvar[n*i+j] = model.getVarByName("mk"+to_string(i)+"_"+to_string(j));
		}
	}

	//Fix the difference
	for(uint i = 0; i < 2*n; i++){
		model.addConstr(plaintext[i] == inputDiff[i]);
		model.addConstr(ciphertext[i] == outputDiff[i]);
	}
	for(uint i = 0; i < m*n; i++)
		model.addConstr(keyvar[i] == keyDiff[i]);

	//Arbitrary objective to help the solver
	GRBLinExpr obj = 0;
	for(uint i = 0; i < 32; i++){
		obj += model.getVarByName("x"+to_string(nbRound/2)+"_"+to_string(i));
		obj += model.getVarByName("y"+to_string(nbRound/2)+"_"+to_string(i));
	}

	model.setObjective(obj, GRB_MINIMIZE);
	model.set(GRB_IntParam_SolutionLimit, 1);
	//Focus on feasibility
	model.set(GRB_IntParam_MIPFocus,1);

	model.update();
	model.optimize();
	return model.get(GRB_IntAttr_SolCount) > 0;
}

void checkTruncatedDiffSpeckRelatedKey(uint const n,
									   uint const m,
									   uint const nbRound,
									   uint const gamma,
									   std::vector<uint> const & inputIndex,
									   std::vector<uint> const & outputIndex,
									   std::vector<uint> const & keyIndex){

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);

	uint64_t boundInput = (1ULL << inputIndex.size());
	uint64_t boundOutput = (1ULL << outputIndex.size());

	uint64_t ctrPossible = 0;
	uint64_t ctrImpossible = 0;

	cout << "---------------------------------------" << endl;
	cout << "Searching for truncated RK RX-diff for Speck with" << endl;
	cout << "n = " << n << " m = " << m << " nbRound = " << nbRound << endl;
	cout << "Input Index = {";
	for(auto const tmp : inputIndex)
		cout << tmp << ",";
	cout << "}" << endl << "Output Index = {";
	for(auto const tmp : outputIndex)
		cout << tmp << ",";
	cout << "}" << endl;
	cout << "Key Index = {";
	for(auto const tmp : keyIndex)
		cout << tmp << ",";
	cout << "}" << endl;
	cout << "---------------------------------------" << endl;

	//Create the key diff
	vector<uint> keyDiff(n*m,0);
	for(auto const tmp : keyIndex)
		keyDiff[tmp] = 1;

	for(uint64_t xinput = 0; xinput < boundInput; xinput++){

		//Create the input diff 
		vector<uint> inputDiff(64,0);
		vector<uint> currentInput;
		for(uint i = 0; i < inputIndex.size(); i++){
			if((xinput >> i)&1){
				inputDiff[inputIndex[i]] = 1;
				currentInput.emplace_back(inputIndex[i]);
			}
		}

		for(uint64_t xoutput = 0; xoutput < boundOutput; xoutput++){

			//Create the output diff
			vector<uint> outputDiff(64,0);
			vector<uint> currentOutput;
			for(uint i = 0; i < outputIndex.size(); i++){
				if((xoutput >> i)&1){
					outputDiff[outputIndex[i]] = 1;
					currentOutput.emplace_back(outputIndex[i]);
				}
			}
				

			bool check = existTrailSpeckRelatedKey(n,m,nbRound,gamma,inputDiff,outputDiff,keyDiff,env);
			
			if(check){
				cout << "For inputDiff = {";
				for(auto const tmp : currentInput)
					cout << tmp << ",";
				cout << "} outputDiff = {";
				for(auto const tmp : currentOutput)
					cout << tmp << ",";
				cout << "} --> ";
				cout << "Trail... :(" << endl;
				ctrPossible++;
			}
			else{
				// cout << "For inputDiff = {";
				// for(auto const tmp : currentInput)
				// 	cout << tmp << ",";
				// cout << "} outputDiff = {";
				// for(auto const tmp : currentOutput)
				// 	cout << tmp << ",";
				// cout << "} --> ";
				// cout << "Impossible Diff !!!" << endl;
				ctrImpossible++;
			}
		}
	}

	cout << ctrPossible << " differentials with a trail" << endl;
	cout << ctrImpossible << " impossible differentials" << endl;
}

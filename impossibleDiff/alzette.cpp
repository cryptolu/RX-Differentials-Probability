#include "alzette.hpp"
using namespace std;

GRBModel getAlzetteModel(std::vector<uint> const & params,
						 uint64_t const cst,
						 uint const gamma,
						 GRBEnv & env){

	uint nbRound = params.size()/2;
	//Write the constant in binary
	vector<uint> bincst(32);
	for(uint i = 0; i < 32; i++)
		bincst[i] = (cst >> i)&1;

	// Create an empty model
    GRBModel model = GRBModel(env);

    //State variables
    vector<vector<GRBVar>> x(nbRound+1, vector<GRBVar>(32));
    vector<vector<GRBVar>> y(nbRound+1, vector<GRBVar>(32));
    vector<vector<GRBVar>> z(nbRound, vector<GRBVar>(32));
    for(uint i = 0; i < nbRound+1; i++){
    	for(uint j = 0; j < 32; j++){
    		x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i)+"_"+to_string(j));
    		y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i)+"_"+to_string(j));
    	}
    }
    for(uint i = 0; i < nbRound; i++){
    	for(uint j = 0; j < 32; j++)
    		z[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z"+to_string(i)+"_"+to_string(j));
    }
    model.update();

    for(uint r = 0;  r < nbRound; r++){

    	vector<GRBVar> tmp(32);

    	//z[r] = x[r] + (y[r] >>> params[2r])
    	for(uint i = 0; i < 32; i++)
    		tmp[i] = y[r][(i+params[2*r])%32]; //tmp = y[r] >>> params[2r]
    	addModAddRXDiffConstr(model, x[r], tmp, z[r], gamma);

    	//x[r+1] = z[r] ^ cst
    	addRXCstXORConstr(model, z[r], bincst, x[r+1], gamma);

    	//y[r+1] = y[r] ^ (z[r] >>> params[2r+1])
    	for(uint i = 0; i < 32; i++)
    		tmp[i] = z[r][(i+params[2*r+1])%32]; //tmp = z[r] >>> params[2r+1]
    	for(uint i = 0; i < 32; i++)
    		addXORConstr(model, y[r][i], tmp[i], y[r+1][i]);
    }

    return model;
}

void searchWeight1RXDiffAlzette(std::vector<uint> const & params,
								uint64_t const cst){

	auto time_start = std::chrono::high_resolution_clock::now();

	cout << "---------------------------------------" << endl;
	cout << "Alzette Search for Weight 1 RXDiff with" << endl;
	cout << "params = {";
	for(auto const tmp : params)
		cout << tmp << ",";
	cout << "}; cst = " << setfill('0') << setw(8) << hex << cst << dec << endl;
	cout << "---------------------------------------" << endl;

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);
	uint n = 32;
	uint nbRound = params.size()/2;
	
	#pragma omp parallel for
	for(uint gamma = 0; gamma < n; gamma++){
		// cout << "------ For gamma = " << gamma << " ------" << endl;
		uint nbPossible = 0;
		uint nbImpossible = 0;
		vector<vector<uint>> listImpossible(2*n);

		for(uint indexInput = 0; indexInput < 2*n; indexInput++){
			for(uint indexOutput = 0; indexOutput < 2*n; indexOutput++){

				// cout << indexInput << "/64 " << indexOutput << "/64            \r" << flush;

				auto model = getAlzetteModel(params,cst,gamma,env);
				vector<GRBVar> plaintext(64);
				vector<GRBVar> ciphertext(64);
				for(uint i = 0; i < 32; i++){
					plaintext[i] = model.getVarByName("x0_"+to_string(i));
					ciphertext[i] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i));
				}
				for(uint i = 32; i < 64; i++){
					plaintext[i] = model.getVarByName("y0_"+to_string(i-32));
					ciphertext[i] = model.getVarByName("y"+to_string(nbRound)+"_"+to_string(i-32));
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
				if(model.get(GRB_IntAttr_SolCount) > 0){
					nbPossible++;
					// for(uint r = 0; r < nbRound+1; r++){
					// 	if(r < 10) cout << " ";
					// 	cout << "x" << r << ": ";
					// 	uint32_t valxi = 0;
					// 	for(uint i = 0; i < 32; i++){
					// 		GRBVar v = model.getVarByName("x"+to_string(r)+"_"+to_string(i));
					// 		valxi |= uint32_t(round(v.get(GRB_DoubleAttr_X))) << i;
					// 	}
					// 	cout << setfill('0') << setw(8) << hex << valxi << " " << dec;

					// 	cout << "y" << r << ": ";
					// 	uint32_t valyi = 0;
					// 	for(uint i = 0; i < 32; i++){
					// 		GRBVar v = model.getVarByName("y"+to_string(r)+"_"+to_string(i));
					// 		valyi |= uint32_t(round(v.get(GRB_DoubleAttr_X))) << i;
					// 	}
					// 	cout << setfill('0') << setw(8) << hex << valyi << " " << dec;
					// 	cout << endl;
					// }
					// cout << "----------------" << endl;
				}
				else{
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

			}
			cout << indexInput << " -> {";
			for(auto const tmp : listImpossible[indexInput])
				cout << tmp << ",";
			cout << "} (" << listImpossible[indexInput].size() << " imp)" << endl;
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

bool existTrailAlzette(std::vector<uint> const & params,
					   uint64_t const cst,
					   uint const gamma,
					   std::vector<uint> const & inputDiff,
					   std::vector<uint> const & outputDiff,
					   GRBEnv & env){

	uint n = 32;
	uint nbRound = params.size()/2;
	auto model = getAlzetteModel(params,cst,gamma,env);
	vector<GRBVar> plaintext(64);
	vector<GRBVar> ciphertext(64);
	for(uint i = 0; i < 32; i++){
		plaintext[i] = model.getVarByName("x0_"+to_string(i));
		ciphertext[i] = model.getVarByName("x"+to_string(nbRound)+"_"+to_string(i));
	}
	for(uint i = 32; i < 64; i++){
		plaintext[i] = model.getVarByName("y0_"+to_string(i-32));
		ciphertext[i] = model.getVarByName("y"+to_string(nbRound)+"_"+to_string(i-32));
	}

	for(uint i = 0; i < 2*n; i++){
		model.addConstr(plaintext[i] == inputDiff[i]);
		model.addConstr(ciphertext[i] == outputDiff[i]);
	}
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

void checkTruncatedDiffAlzette(std::vector<uint> const & params,
							   uint64_t const cst,
							   uint const gamma,
							   vector<uint> const & inputIndex,
							   vector<uint> const & outputIndex){

	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
	env.set(GRB_IntParam_Threads, 1);

	uint64_t boundInput = (1ULL << inputIndex.size());
	uint64_t boundOutput = (1ULL << outputIndex.size());

	uint64_t ctrPossible = 0;
	uint64_t ctrImpossible = 0;

	cout << "---------------------------------------" << endl;
	cout << "Searching for truncated RX-diff for Alzette with" << endl;
	cout << "params = {";
	for(auto const tmp : params)
		cout << tmp << ",";
	cout << "}; cst = " << setfill('0') << setw(8) << hex << cst << dec << endl;
	cout << "Input Index = {";
	for(auto const tmp : inputIndex)
		cout << tmp << ",";
	cout << "}" << endl << "Output Index = {";
	for(auto const tmp : outputIndex)
		cout << tmp << ",";
	cout << "}" << endl;
	cout << "gamma = " << gamma << endl;
	cout << "---------------------------------------" << endl;

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

			bool check = existTrailAlzette(params,cst,gamma,inputDiff,outputDiff,env);
			
			if(check){
				// cout << "For inputDiff = {";
				// for(auto const tmp : currentInput)
				// 	cout << tmp << ",";
				// cout << "} outputDiff = {";
				// for(auto const tmp : currentOutput)
				// 	cout << tmp << ",";
				// cout << "} --> ";
				// cout << "Trail... :(" << endl;
				ctrPossible++;
			}
			else{
				cout << "For inputDiff = {";
				for(auto const tmp : currentInput)
					cout << tmp << ",";
				cout << "} outputDiff = {";
				for(auto const tmp : currentOutput)
					cout << tmp << ",";
				cout << "} --> ";
				cout << "Impossible Diff !!!" << endl;
				ctrImpossible++;
			}
		}
	}

	cout << ctrPossible << " differentials with a trail" << endl;
	cout << ctrImpossible << " impossible differentials" << endl;

}
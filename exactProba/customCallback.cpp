#include "customCallback.hpp"

using namespace std;

typedef unsigned int uint;

CustomCallback::CustomCallback(uint const wordSize,
							   uint const rotationOffset,
							   std::vector<std::array<std::vector<GRBVar>,3>> const & modAddVars,
							   std::vector<std::vector<GRBVar>> const & stateVars,
							   GRBVar & sumneq,
							   std::vector<std::vector<GRBVar>> const & keyVars):
							   wordSize(wordSize),
							   rotationOffset(rotationOffset),
							   modAddVars(modAddVars),
							   stateVars(stateVars),
							   sumneq(sumneq),
							   keyVars(keyVars)
{}

void CustomCallback::callback(){
try{
	if (where == GRB_CB_MIPSOL) { //If the solver found a solution
		
		cout << "in callback " << endl;
		uint nbModAdd = modAddVars.size();
		cout << "nbModAdd = " << nbModAdd << endl;
		//Get the value of the mod add alpha,beta,delta
		vector<array<uint64_t,3>> val_abd(nbModAdd);
		for(uint r = 0; r < nbModAdd; r++){
			for(uint i = 0; i < 3; i++){
				uint64_t wordVal = 0;
				for(uint j = 0; j < wordSize; j++){
					uint tmp = uint(round(getSolution(modAddVars[r][i][j])));
					if(tmp)
						wordVal |= (1ULL << j);
				}
				val_abd[r][i] = wordVal;
			}
		}


		//Get the value of the state for the trail
		uint nbRound = stateVars.size();
		cout << "nbRound = " << nbRound << endl;
		vector<vector<uint64_t>> val_state(nbRound);
		for(uint r = 0; r < nbRound; r++){
			uint nbWord = stateVars[r].size()/wordSize;
			for(uint i = 0; i < nbWord; i++){
				uint64_t wordVal = 0;
				for(uint64_t j = 0; j < wordSize; j++){
					uint tmp = uint(round(getSolution(stateVars[r][i*wordSize + j])));
					if(tmp)
						wordVal |= (1ULL << j);
				}
				val_state[r].emplace_back(wordVal);
			}
		}

		//Get the value for the keyvars, if any
		uint nbKey = keyVars.size();
		vector<vector<uint64_t>> val_key(nbKey);
		for(uint r = 0;  r < nbKey; r++){
			uint nbWord = keyVars[r].size()/wordSize;
			for(uint i = 0; i < nbWord; i++){
				uint64_t wordVal = 0;
				for(uint64_t j = 0; j < wordSize; j++){
					uint tmp = uint(round(getSolution(keyVars[r][i*wordSize + j])));
					if(tmp)
						wordVal |= (1ULL << j);
				}
				val_key[r].emplace_back(wordVal);
			}
		}

		uint val_sumneq = uint(round(getSolution(sumneq)));

		//Print the trail
		cout << "=========================================================" << endl;
		cout << "======== Solution found with Objective : " << getDoubleInfo(GRB_CB_MIPSOL_OBJ) << " ========" << endl;
		cout << "Number of notAllEqual vars : " << val_sumneq << endl;
		cout << "Trail :" << endl;
		for(uint r = 0; r < nbRound; r++){
			cout << "Round " << r << " : ";
			for(auto const & word : val_state[r])
				cout << "0x" << setfill('0') << setw(wordSize/4) << hex << word << dec << " ";
			if(r < val_key.size()){
				cout << " k" << r << " ";
				for(auto const & word : val_key[r])
					cout << "0x" << setfill('0') << setw(wordSize/4) << hex << word << dec << " ";
			}
			cout << endl;
		}

		//Print the mod add transitions and associated probabilities
		double sumLog = 0;
		for(uint r = 0; r < nbModAdd; r++){
			//Count the number of solutions, deduce the log
			uint64_t countsol = getRXDiffCount(val_abd[r][0],val_abd[r][1],val_abd[r][2],wordSize,rotationOffset);
			auto valLog = log2(countsol)-2*wordSize;
			sumLog += valLog;
			cout << "0x" << setfill('0') << setw(wordSize/4) << hex << val_abd[r][0];
			cout << " + ";
			cout << "0x" << setfill('0') << setw(wordSize/4) << hex << val_abd[r][1];
			cout << " -> ";
			cout << "0x" << setfill('0') << setw(wordSize/4) << hex << val_abd[r][2];
			cout << dec;
			cout << " (2^" << valLog << ")" << endl;
		}
		cout << "sumLog = " << sumLog << endl;

		//Print the modAddTrail for easy sage/python copypaste
		cout << "modAddTrail = [";
		for(uint r = 0; r < nbModAdd; r++){
			cout << "[";
			cout << "0x" << setfill('0') << setw(wordSize/4) << hex << val_abd[r][0];
			cout << ",";
			cout << "0x" << setfill('0') << setw(wordSize/4) << hex << val_abd[r][1];
			cout << ",";
			cout << "0x" << setfill('0') << setw(wordSize/4) << hex << val_abd[r][2];
			cout << dec;
			if(r < nbModAdd-1)
				cout << "]," << endl;
			else
				cout << "]]" << endl;
		}
		cout << "=========================================================" << endl;
	}
} catch (GRBException e) {
cout << "Error number: " << e.getErrorCode() << endl;
cout << e.getMessage() << endl;
} catch (...) {
cout << "Error during callback" << endl;
}
}






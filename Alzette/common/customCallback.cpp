#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "customCallback.hpp"
#include "../common/rxdp.hpp"

using namespace std;

typedef unsigned int uint;

char *filename = NULL;
string current_cipher = "unknown";

void setCipherName(char * cipher) {
	assert(strlen(cipher) > 0);
	current_cipher = cipher;
}

void setFilename(char * _filename) {
	filename = _filename;
	if (filename == NULL) {
		return;
	}
	FILE * fd;
	fd = fopen(filename, "r");
	if (fd) {
		printf("file %s already exists! remove it to start from scratch\n", filename);
		exit(-1);
	}
	fd = fopen(filename, "w");
	assert(fd);
	fclose(fd);
}

// struct Trail {
// string cipher;
//	uint wordSize;
//	uint rotationOffset;

// uint nbStateWords;
// uint nbKeyWords;
	
// 	uint sumNeq;
// 	double sumLog;

// 	std::vector<std::vector<uint64_t>> stateWords;
// 	std::vector<std::vector<uint64_t>> keyWords;
// 	std::vector<std::vector<uint64_t>> modAddWords;

// 	void write_to_file(char *filename);
// };
void Trail::write_to_file(char *filename) {
	FILE *fd = fopen(filename, "w");
	assert(fd);

	fprintf(fd, "%s\n", cipher.c_str());

	fprintf(fd, "%u %u %u %u\n", wordSize, rotationOffset, nbStateWords, nbKeyWords);
	fprintf(fd, "%u %lf\n", sumNeq, sumLog);

	fprintf(fd, "%lu\n", stateWords.size());
	for (auto &vec: stateWords) {
		assert(vec.size() == nbStateWords);
		for (auto &word: vec)
			fprintf(fd, "%08lx ", word);
		fprintf(fd, "\n");
	}

	fprintf(fd, "%lu\n", keyWords.size());
	for (auto &vec: keyWords) {
		assert(vec.size() == nbKeyWords);
		for (auto &word: vec)
			fprintf(fd, "%08lx ", word);
		fprintf(fd, "\n");
	}

	fprintf(fd, "%lu\n", modAddWords.size());
	for (auto &vec: modAddWords) {
		assert(vec.size() == 3);
		for (auto &word: vec)
			fprintf(fd, "%08lx ", word);
		fprintf(fd, "\n");
	}
	
	fclose(fd);
}

void Trail::read_from_file(char *filename) {
	size_t vec_size;

	FILE *fd = fopen(filename, "r");
	assert(fd);

	{
		char buf[4096] = {};
		assert(1 == fscanf(fd, "%s\n", buf));
		cipher = string(buf);
	}

	assert(4 == fscanf(fd, "%u %u %u %u\n", &wordSize, &rotationOffset, &nbStateWords, &nbKeyWords));
	assert(2 == fscanf(fd, "%u %lf\n", &sumNeq, &sumLog));
	if (!wordSize) {
		printf("Trail file %s corrupted!\n", filename);
		exit(-1);
	}

	stateWords.clear();
	assert(1 == fscanf(fd, "%lu\n", &vec_size));
	for (size_t i = 0; i < vec_size; i++) {
		vector<uint64_t> vec;
		for (size_t j = 0; j < nbStateWords; j++) {
			uint64_t word;
			assert(1 == fscanf(fd, "%08lx ", &word));
			vec.push_back(word);
		}
		stateWords.push_back(vec);
		fscanf(fd, "\n");
	}

	keyWords.clear();
	assert(1 == fscanf(fd, "%lu\n", &vec_size));
	for (size_t i = 0; i < vec_size; i++) {
		vector<uint64_t> vec;
		for (size_t j = 0; j < nbKeyWords; j++) {
			uint64_t word;
			assert(1 == fscanf(fd, "%08lx ", &word));
			vec.push_back(word);
		}
		keyWords.push_back(vec);
		fscanf(fd, "\n");
	}

	modAddWords.clear();
	assert(1 == fscanf(fd, "%lu\n", &vec_size));
	for (size_t i = 0; i < vec_size; i++) {
		vector<uint64_t> vec;
		for (size_t j = 0; j < 3; j++) {
			uint64_t word;
			assert(1 == fscanf(fd, "%08lx ", &word));
			vec.push_back(word);
		}
		modAddWords.push_back(vec);
		fscanf(fd, "\n");
	}

	fclose(fd);
}

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
		Trail trail;

		trail.cipher = current_cipher;

		trail.wordSize = wordSize;
		trail.rotationOffset = rotationOffset;

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
		assert(nbRound > 0);
		trail.nbStateWords = val_state[0].size();
		trail.stateWords = val_state;


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
		if (nbKey > 0) {
			assert(val_key.size());
			trail.nbKeyWords = val_key[0].size();
		}
		else {
			trail.nbKeyWords = 0;
		}
		trail.keyWords = val_key;

		uint val_sumneq = uint(round(getSolution(sumneq)));
		trail.sumNeq = val_sumneq;

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
		vector<vector<uint64_t>> val_modAdd;
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
			val_modAdd.push_back({val_abd[r][0], val_abd[r][1], val_abd[r][2]});
		}
		cout << "sumLog = " << sumLog << endl;
		trail.sumLog = sumLog;
		trail.modAddWords = val_modAdd;

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

		if (filename)
			trail.write_to_file(filename);
	}
} catch (GRBException e) {
cout << "Error number: " << e.getErrorCode() << endl;
cout << e.getMessage() << endl;
} catch (...) {
cout << "Error during callback" << endl;
}
}






#include "speck.hpp"
using namespace std;

std::array<std::vector<std::vector<GRBVar>>,2> 
addSpeckKSConstr(GRBModel & model,
				 uint const n,
				 uint const m,
				 uint const nbRound,
				 uint const alpha,
				 uint const beta,
				 std::string const & suffix){
//Add constraints for the key schedule of Speck
//These are constraints *in values* not in differences, so constants have an impact
//Return the master key and round key variables

	//Create variables for the key schedule
	auto mk = genArray2DBinVar(model, m, n, "mk"+suffix);
    auto k = genArray2DBinVar(model, nbRound, n, "k"+suffix);
	//We need the l[i] variables from the specification, but also some additional temp variables
	auto l = genArray2DBinVar(model, nbRound+m-2, n, "l"+suffix);
	auto tmpKS = genArray2DBinVar(model, nbRound-1, n, "tmpKS"+suffix);
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
	for(uint64_t i = 0; i < nbRound-1; i++){
		vector<GRBVar> Sali(n); //l[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Sali[j] = l[i][mod(j+alpha,n)];

		vector<GRBVar> Sbki(n); //k[i] rotated to the left by beta
		for(uint j = 0; j < n; j++)
			Sbki[j] = k[i][mod(j-beta,n)];

		//l[i+m-1] = (k[i] + Sali) ^ i
		//First, tmp = k[i] + Sali
		addModAddValueConstr(model, k[i], Sali, tmpKS[i], "carryKS"+suffix+to_string(i));

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

	array<vector<vector<GRBVar>>,2> ret = {mk,k};
    return ret;
}

std::array<std::vector<std::vector<GRBVar>>,3> 
addSpeckRXDConstr(GRBModel & model,
				  uint const n,
				  uint const nbRound,
				  uint const alpha, 
				  uint const beta, 
				  uint const gamma,
				  std::vector<std::vector<GRBVar>> & k){
//Add constraints for RX-differential propagation in Speck
//Returns a triplet for containing the variables [x,y,z]

	//State variables
	auto x = genArray2DBinVar(model, nbRound+1, n, "dx");
	auto y = genArray2DBinVar(model, nbRound+1, n, "dy");
	auto z = genArray2DBinVar(model, nbRound, n, "dz");

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

		//y[i+1] = Sbyi ^ x[i+1]
		for(uint j = 0; j < n; j++)
			addXORConstr(model, Sbyi[j], x[i+1][j], y[i+1][j]);

		//xor with the key
    	//x[i+1] = z[i] ^ k[i]
		addRXCstXORConstr(model, z[i], k[i], x[i+1], gamma);
	}

	array<vector<vector<GRBVar>>,3> ret = {x,y,z};
    return ret;
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
//Returns the model

	// Create an empty model
    GRBModel model = GRBModel(env);

    //Key schedule constraints
	auto [mk,k] = addSpeckKSConstr(model,n,m,nbRound,alpha,beta);

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

	//Constraints for the differential trail
    auto [x,y,z] = addSpeckRXDConstr(model,n,nbRound,alpha,beta,gamma,k);
	

	//notAllEqual variables
    vector<vector<GRBVar>> v(nbRound, vector<GRBVar>(n-1));
    for(uint i = 0; i < nbRound; i++){
    	for(uint j = 0; j < gamma-1; j++)
    		v[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "neq"+to_string(i)+"_"+to_string(j));
    	//don't check the bit of index gamma-1, as it's the MSB of the right half
    	for(uint j = gamma; j < n-1; j++)
    		v[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "neq"+to_string(i)+"_"+to_string(j));
    }
    //Sum of neq vars in the same round
    vector<GRBVar> sumneq_round(nbRound);
    for(uint r = 0; r < nbRound; r++)
		sumneq_round[r] = model.addVar(0.0, n, 0.0, GRB_INTEGER, "sumneq_round"+to_string(r));
    //Sum of all neq vars
    GRBVar sumneq = model.addVar(0.0, nbRound*n, 0.0, GRB_INTEGER, "sumneq");
    model.update();


    GRBLinExpr expr_sumneq = 0;
	for(uint r = 0; r < nbRound; r++){
		GRBLinExpr expr_sumneq_round = 0;
		for(uint i = 0; i < gamma-1; i++)
			expr_sumneq_round += v[r][i];
		for(uint i = gamma; i < n-1; i++)
			expr_sumneq_round += v[r][i];
		model.addConstr(expr_sumneq_round == sumneq_round[r]);
		expr_sumneq += sumneq_round[r];
	}
	model.addConstr(expr_sumneq == sumneq);

	for(uint i = 0; i < nbRound; i++){
		vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Saxi[j] = x[i][mod(j+alpha,n)];

		//notAllEq constraints
    	for(uint j = 0; j < gamma-1; j++){
    		vector<GRBVar> vars({Saxi[j], y[i][j], z[i][j]});
    		addNotAllEqualBoolConstr(model, vars, v[i][j]);
    	}
    	for(uint j = gamma; j < n-1; j++){
    		vector<GRBVar> vars({Saxi[j], y[i][j], z[i][j]});
    		addNotAllEqualBoolConstr(model, vars, v[i][j]);
    	}
	}

	return model;
}

std::vector<GRBVar> 
addModAddProbaSpeck(GRBModel & model,
					uint const alpha,
					uint const gamma,
					std::vector<std::vector<GRBVar>> & x,
					std::vector<std::vector<GRBVar>> & y,
					std::vector<std::vector<GRBVar>> & z){
//Add mod add probability constraints and variables
//Return a vector containing the weight variables
	
	uint nbRound = x.size()-1;
	uint n = x[0].size();

	//Weight variables
    vector<GRBVar> wt(nbRound);
    for(uint i = 0; i < nbRound; i++){
    	wt[i] = model.addVar(0,2*n,0,GRB_CONTINUOUS, "wt"+to_string(i));
    }
    model.update();
	//Round functions
	for(uint i = 0; i < nbRound; i++){
		vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Saxi[j] = x[i][mod(j+alpha,n)];

    	//Probability constraints
    	addModAddRXProbaConstr(model,Saxi,y[i],z[i],wt[i],gamma,"_r"+to_string(i));
	}
	model.update();

	return wt;
}

std::tuple<double,uint, std::vector<std::array<uint64_t,3>>> 
searchBestTrailSpeck(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				bool const withExactProbability,
				std::map<uint,uint> const & lowerBoundConsecutiveRounds,
				std::vector<uint64_t> const & inputDiff,
				std::vector<uint64_t> const & outputDiff){

	uint val_sumneq = 0;
	vector<array<uint64_t,3>> trail;
	auto knownTrail = make_pair(val_sumneq, trail);
	return searchBestTrailSpeck(params,gamma,env,withExactProbability,lowerBoundConsecutiveRounds,knownTrail,inputDiff,outputDiff);
}

std::tuple<double,uint, std::vector<std::array<uint64_t,3>>> 
searchBestTrailSpeck(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				bool const withExactProbability,
				std::map<uint,uint> const & lowerBoundConsecutiveRounds,
				std::pair<uint, std::vector<std::array<uint64_t,3>>> const & knownTrail,
				std::vector<uint64_t> const & inputDiff,
				std::vector<uint64_t> const & outputDiff){
/*
Search for the best RX-diff trail with RX rotation gamma
params should be the Speck parameters, in order: n,m,nbRound,alpha,beta
knownTrail can be provided as a pair (W,trail) where W is a bound on the number of neq variables and trail is a starting solution for the model
withExactProbability allows to control if we add probability constraints and optimization
False means we just optimize on the number of neq variables
lowerBoundConsecutiveRounds provides lower bounds on the number of Neq vars for consecutive rounds
e.g. if lowerBoundConsecutiveRounds[2] = 10, then we add constraints so that for any 2 consecutive rounds, there are at least 10 neq vars set to 1
inputDiff and outputDiff can be left empty to allow any input/output differential, or given to fix the input/output differential respectively, as {x,y} for the difference on the left,right branch. Won't work for n > 64 since differences are given as a uint64_t for each word (but shouldn't really be a problem).
knownTrail can be used to help the solver with a starting point. The first element of the pair is an upper bound on the number of neq vars, and the second element of the pair is the trail given as triplets {x,y,k} for the difference on x,y and the value of k at a given round
*/

	auto [n,m,nbRound,alphaRot,betaRot] = params;
	auto model = getSpeckModel(n,m,nbRound,alphaRot,betaRot,gamma,env);

	//State variables
	auto x = get2DArrayVar(model, nbRound+1, n, "dx");
    auto y = get2DArrayVar(model, nbRound+1, n, "dy");
    auto z = get2DArrayVar(model, nbRound, n, "dz");
    auto mk = get2DArrayVar(model, m, n, "mk");
    auto k = get2DArrayVar(model, nbRound, n, "k");

	//Input and output constraints (if any)
	if(inputDiff.size() > 0){
		addValueToBinVectorConstr(model, x[0], inputDiff[0]);
		addValueToBinVectorConstr(model, y[0], inputDiff[1]);
	}
	if(outputDiff.size() > 0){
		addValueToBinVectorConstr(model, x[nbRound], outputDiff[0]);
		addValueToBinVectorConstr(model, y[nbRound], outputDiff[1]);
	}

    //Sum of all neq vars
    GRBVar sumneq = model.getVarByName("sumneq");

	//Prepare the callback
	vector<array<vector<GRBVar>,3>> modAddVars(nbRound);
	for(uint r = 0; r < nbRound; r++){
		modAddVars[r][0] = vector<GRBVar>(n);
		modAddVars[r][1] = vector<GRBVar>(n);
		modAddVars[r][2] = vector<GRBVar>(n);
		for(uint i = 0; i < n; i++){
			modAddVars[r][0][i] = x[r][mod(i+alphaRot,n)];
			modAddVars[r][1][i] = y[r][i];
			modAddVars[r][2][i] = z[r][i];
		}
	}
	vector<vector<GRBVar>> stateVars(nbRound+1, vector<GRBVar>(2*n));
	for(uint r = 0; r < nbRound+1; r++){
		for(uint i = 0; i < n; i++){
			stateVars[r][i] = x[r][i];
			stateVars[r][i+n] = y[r][i];
		}
	}
	model.update();
	CustomCallback cb(n,gamma,modAddVars,stateVars,sumneq,k);
	model.setCallback(&cb);

	//Check if we have a valid trail provided as input, if so, use it to bound the number of neq variables and provide a starting value for the trail variables
	auto const & trail = knownTrail.second;
    if(trail.size() > 0){
		//Bound the number of notAllEq
		// auto const bestBound = knownTrail.first;
		// model.addConstr(sumneq <= bestBound);

		//Provide MIP Start info from the trail
		for(uint r = 0; r < nbRound+1; r++){
			auto const & trail_r = trail[r];
			for(uint i = 0; i < n; i++){
				uint val = (trail_r[0] >> i)&1;
				x[r][i].set(GRB_DoubleAttr_Start, val);

				val = (trail_r[1] >> i)&1;
				y[r][i].set(GRB_DoubleAttr_Start, val);

				if(r < nbRound){
					val = (trail_r[2] >> i)&1;
					k[r][i].set(GRB_DoubleAttr_Start, val);
				}
			}
		}
    }

    //Lower bounds based on results on a smaller number of rounds
    if(lowerBoundConsecutiveRounds.size() > 0){
    	//Get the variable representing the sum of neq vars over a given round
    	auto sumneq_round = getArrayVar(model,nbRound,"sumneq_round");
    	addMatsuiLikeConstr(model, lowerBoundConsecutiveRounds, sumneq_round, sumneq);
    }

    //Set the objective depending on withExactProbability
    if(withExactProbability){
    	auto wt = addModAddProbaSpeck(model,alphaRot,gamma,x,y,z);
    	GRBLinExpr obj = 0;
		for(uint i = 0; i < nbRound; i++)
			obj += wt[i];
		//Objective: minimize probability
		model.setObjective(obj, GRB_MINIMIZE);
    }
    else{
    	//Objective: minimize the number of notAllEq
		GRBLinExpr obj = sumneq;
		model.setObjective(obj, GRB_MINIMIZE);
    }

	model.update();
	model.optimize();
	if(model.get(GRB_IntAttr_SolCount) > 0){

	    vector<uint64_t> valx(nbRound+1,0);
	    vector<uint64_t> valy(nbRound+1,0);
	    vector<uint64_t> valz(nbRound,0);
	    vector<uint64_t> valk(nbRound,0);
	    vector<uint64_t> valmk(m,0);
	    vector<uint64_t> alpha(nbRound,0);
	    vector<uint64_t> beta(nbRound,0);
	    vector<uint64_t> delta(nbRound,0);
	    vector<double> modAddLog(nbRound,0);
	    double sumLog = 0;
	    for(uint i = 0; i < nbRound; i++){

	    	//State values
	    	valx[i] = getUintSolutionFromBinVector(x[i]);
	    	valy[i] = getUintSolutionFromBinVector(y[i]);
	    	valz[i] = getUintSolutionFromBinVector(z[i]);
	    	valk[i] = getUintSolutionFromBinVector(k[i]);

	    	//zi = (xi >>> alphaRot) + yi
	    	vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
			for(uint j = 0; j < n; j++)
				Saxi[j] = x[i][mod(j+alphaRot,n)];
			uint64_t tmpint = getUintSolutionFromBinVector(Saxi);

			//Keep the differences used in modular additions
	    	alpha[i] = tmpint;
	    	beta[i] = valy[i];
	    	delta[i] = valz[i];

	    	uint64_t countsol = getRXDiffCount(alpha[i],beta[i],delta[i],n,gamma);
			sumLog += (log2(countsol)-2*n);
			modAddLog[i] = log2(countsol)-2*n;
		}

		valx[nbRound] = getUintSolutionFromBinVector(x[nbRound]);
	    valy[nbRound] = getUintSolutionFromBinVector(y[nbRound]);

	    for(uint i = 0; i < m; i++)
			valmk[i] = getUintSolutionFromBinVector(mk[i]);

		uint val_sumneq = uint(round(sumneq.get(GRB_DoubleAttr_X)));
		auto bestObj = model.get(GRB_DoubleAttr_ObjVal);
		auto bestLB = model.get(GRB_DoubleAttr_ObjBound);

		cout << "***********************************************************************" << endl;
		cout << "*** Speck " << 2*n << "/" << m*n << " " << nbRound << " rounds with k = " << gamma << " ***" << endl;
		cout << "*** Best solution found with objective " << bestObj << " ***" << endl;
		if(model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
			cout << "Solution obtained after time-out" << endl;
		else
			cout << "Solution obtained after full optimization" << endl;

		cout << "Number of notAllEqual vars : " << val_sumneq << endl;
		cout << "Best lower bound on objective : " << bestLB << endl;

		cout << "Trail :" << endl;
		cout << "Master key : ";
		for(uint i = 0; i < m; i++)
			cout << "0x" << setfill('0') << setw(n/4) << hex << valmk[i] << " ";
		cout << endl;
		for(uint i = 0; i < nbRound+1; i++){
			cout <<  "x" << i << " : " << "0x" << setfill('0') << setw(n/4) << hex << valx[i];
			cout << " y" << i << " : " << "0x" << setfill('0') << setw(n/4) << hex << valy[i];
			if(i < nbRound)
				cout << " k" << i << " : " << "0x" << setfill('0') << setw(n/4) << hex << valk[i];
			cout << endl;
		}
		cout << "modAdd :" << endl;
		for(uint i = 0; i < nbRound; i++){
			cout << "0x" << setfill('0') << setw(n/4) << hex << alpha[i] << " + ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << beta[i] << " -> ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << delta[i] << dec;
			cout << " (2^" << modAddLog[i] << ")" << endl;
		}
		cout << "sumLog = " << sumLog << endl;
		cout << "modAddTrail = [";
		for(uint i = 0; i < nbRound; i++){
			cout << "[";
			cout << "0x" << setfill('0') << setw(n/4) << hex << alpha[i];
			cout << ",";
			cout << "0x" << setfill('0') << setw(n/4) << hex << beta[i];
			cout << ",";
			cout << "0x" << setfill('0') << setw(n/4) << hex << delta[i];
			if(i < nbRound-1)
				cout << "]," << endl;
			else
				cout << "]]" << endl;
		}
		cout << "***********************************************************************" << endl;
		cout << dec;


		vector<array<uint64_t,3>> trail(nbRound+1);
		for(uint r = 0; r < nbRound+1; r++){
			trail[r][0] = valx[r];
			trail[r][1] = valy[r];
			if(r < nbRound)
				trail[r][2] = valk[r];
		}

		return make_tuple(bestLB, val_sumneq, trail);
	}
	else{
		cout << "No solution" << endl;
		uint val_sumneq = 0;
		double bestLB = 0;
		vector<array<uint64_t,3>> trail;
		return make_tuple(bestLB, val_sumneq, trail);
	}
}

std::array<std::vector<std::vector<GRBVar>>,3> 
addSpeckPlaintextValueConstr(GRBModel & model,
							 uint const n,
							 uint const nbRound,
							 uint const alpha, 
							 uint const beta, 
							 std::string const & suffix,
							 std::vector<std::vector<GRBVar>> & k){
//Add constraints for one plaintext in value for Speck
// suffix is used to create unique variable names
//Give the key variables as input, but keyschedule should be handled separately
//Returns a triplet for containing the variables [x,y,z]

	auto x = genArray2DBinVar(model, nbRound+1, n, "x"+suffix);
	auto y = genArray2DBinVar(model, nbRound+1, n, "y"+suffix);
	auto z = genArray2DBinVar(model, nbRound, n, "z"+suffix);

	for(uint i = 0; i < nbRound; i++){
		vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Saxi[j] = x[i][mod(j+alpha,n)];

		vector<GRBVar> Sbyi(n); //k[i] rotated to the left by beta
		for(uint j = 0; j < n; j++)
			Sbyi[j] = y[i][mod(j-beta,n)];

		//zi = Saxi + yi
		addModAddValueConstr(model,Saxi,y[i],z[i],"carry"+suffix+to_string(i));

		//x[i+1] = z[i] ^ k[i]
		//y[i+1] = Sbyi ^ x[i+1]
		for(uint j = 0; j < n; j++){
			addXORConstr(model, z[i][j], k[i][j], x[i+1][j]);
			addXORConstr(model, Sbyi[j], x[i+1][j], y[i+1][j]);
		}
	}

	array<vector<vector<GRBVar>>,3> ret = {x,y,z};
    return ret;
}

GRBModel
getSpeckModelPlaintextPair(uint const n,
						   uint const m, 
						   uint const nbRound,
						   uint const alpha, 
						   uint const beta, 
						   uint const gamma,
						   GRBEnv & env,
						   std::vector<std::array<uint64_t,3>> const & trail){
//Create a MILP model for Speck to find a plaintext pair matching a given trail
	// Create an empty model
    GRBModel model = GRBModel(env);

    //Key schedule constraints
	auto [mk,k] = addSpeckKSConstr(model,n,m,nbRound,alpha,beta);

	//Constraints for the differential trail
    auto [x,y,z] = addSpeckRXDConstr(model,n,nbRound,alpha,beta,gamma,k);

    //Constraints for the plaintext values
    auto [x0,y0,z0] = addSpeckPlaintextValueConstr(model,n,nbRound,alpha,beta,"0",k);
    auto [x1,y1,z1] = addSpeckPlaintextValueConstr(model,n,nbRound,alpha,beta,"1",k);

    //--- Bind the value variables to the differential variables ---
    //Essentially, enforce dx = rot(x0) xor x1
    for(uint r = 0; r < nbRound+1; r++){
    	addRXDifferentialValueBinding(model,x0[r],x1[r],x[r],gamma);
    	addRXDifferentialValueBinding(model,y0[r],y1[r],y[r],gamma);
    }
    for(uint r = 0; r < nbRound; r++)
    	addRXDifferentialValueBinding(model,z0[r],z1[r],z[r],gamma);

    //Add the constraint values for the trail
    //Fixing the value of (x,y,k) is enough to deduce a unique value for z
    //The solver should handle that easily during presolving
    for(uint r = 0; r < nbRound+1; r++){
    	addValueToBinVectorConstr(model,x[r],trail[r][0]);
    	addValueToBinVectorConstr(model,y[r],trail[r][1]);
    }
    for(uint r = 0; r < nbRound; r++)
    	addValueToBinVectorConstr(model,k[r],trail[r][2]);

    //Some parameters that should help solving
    //We're only interested in getting one solution
    model.set(GRB_IntParam_SolutionLimit, 1);
    //Focus on finding solutions rather than optimization
    model.set(GRB_IntParam_MIPFocus, 1);
    //Still give an arbitrary objective as it's supposed to help the solver a bit
    GRBLinExpr obj = 0;
    for(uint i = 0; i < n; i++){
    	obj += x0[0][i];
    	obj += y0[0][i];
    }
    model.setObjective(obj, GRB_MINIMIZE);
    return model;
}

std::pair<bool, std::array<uint64_t, 4>>
findPlaintextForTrailSpeck(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				std::vector<std::array<uint64_t,3>> const & trail){
/*
//Return a pair {b, p}, b = true if a plaintext pair was found (and then p contains the pair as x0,y0,x1,y1)
*/

	auto [n,m,nbRound,alpha,beta] = params;
	auto model = getSpeckModelPlaintextPair(n,m,nbRound,alpha,beta,gamma,env,trail);
	model.optimize();

	//Extract the solution, if any
	if(model.get(GRB_IntAttr_SolCount) == 0){
		cout << "*************************************" << endl;
		cout << "*** No matching pair of plaintext ***" << endl;
		cout << "*************************************" << endl;
		array<uint64_t,4> res;
		return make_pair(false,res);
	}
	else{

		auto x0 = getArrayVar(model,n,"x00_");
		auto x1 = getArrayVar(model,n,"x10_");
		auto y0 = getArrayVar(model,n,"y00_");
		auto y1 = getArrayVar(model,n,"y10_");

		array<uint64_t,4> res;
		res[0] = getUintSolutionFromBinVector(x0);
		res[1] = getUintSolutionFromBinVector(y0);
		res[2] = getUintSolutionFromBinVector(x1);
		res[3] = getUintSolutionFromBinVector(y1);

		//Get the round keys from the trail directly
		vector<uint64_t> k(nbRound);
		for(uint r = 0; r < nbRound; r++)
			k[r] = trail[r][2];

		uint64_t dx = CSHL(res[0],gamma,n) ^ res[2];
		uint64_t dy = CSHL(res[1],gamma,n) ^ res[3];
		auto [cx0,cy0] = speckEncrypt(n,nbRound,alpha,beta,res[0],res[1],k);
		auto [cx1,cy1] = speckEncrypt(n,nbRound,alpha,beta,res[2],res[3],k);
		uint64_t cdx = CSHL(cx0,gamma,n) ^ cx1;
		uint64_t cdy = CSHL(cy0,gamma,n) ^ cy1;

		cout << "*******************************" << endl;
		cout << "*** Found pair of plaintext ***" << endl;
		cout << "x0 = 0x" << setfill('0') << setw(n/4) << hex << res[0] << " ";
		cout << "y0 = 0x" << setfill('0') << setw(n/4) << hex << res[1] << endl;
		cout << "x1 = 0x" << setfill('0') << setw(n/4) << hex << res[2] << " ";
		cout << "y1 = 0x" << setfill('0') << setw(n/4) << hex << res[3] << endl;

		//Check the found solution using the actual speck encryption
		//Only check input and output differentials
		cout << "Resulting input dx dy : ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << dx << " ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << dy << endl;
		cout << " Expected input dx dy : ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << trail[0][0] << " ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << trail[0][1] << endl;

		cout << "Resulting output dx dy : ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << cdx << " ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << cdy << endl;
		cout << " Expected output dx dy : ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << trail[nbRound][0] << " ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << trail[nbRound][1] << endl;

		if(dx == trail[0][0] && dy == trail[0][1])
			cout << "Input differentials matches" << endl;
		else
			cout << "Input differentials fail" << endl;
		if(cdx == trail[nbRound][0] && cdy == trail[nbRound][1])
			cout << "Output differentials matches" << endl;
		else
			cout << "Output differentials fail" << endl;

		cout << "*******************************" << endl;

		cout << dec;
		return make_pair(true,res);
	}
}

std::vector<uint64_t>
speckKS(uint const n,
		uint const nbRound,
		uint const alpha,
		uint const beta,
		std::vector<uint64_t> const & mk){

	uint64_t wordMask = (1ULL << n)-1;
	uint m = mk.size();

	vector<uint64_t> k(nbRound);
	vector<uint64_t> l(nbRound+m-2);
	
	k[0] = mk[0];
	for(uint i = 0; i < m-1; i++)
		l[i] = mk[i+1];

	for(uint64_t i = 0; i < nbRound-1; i++){
		uint64_t Sali = CSHR(l[i],alpha,n);
		uint64_t Sbki = CSHL(k[i],beta,n);
		l[i+m-1] = ((k[i] + Sali)&wordMask) ^ i;
		k[i+1] = Sbki ^ l[i+m-1];
	}
	return k;
}

std::array<uint64_t,2>
speckEncrypt(uint const n,
			 uint const nbRound,
			 uint const alpha,
			 uint const beta,
			 uint64_t x,
			 uint64_t y,
			 std::vector<uint64_t> const & k){

	uint64_t wordMask = (1ULL << n)-1;
	for(uint i = 0; i < nbRound; i++){
		x = ((CSHR(x,alpha,n) + y)&wordMask) ^ k[i];
		y = CSHL(y,beta,n) ^ x;
	}
	array<uint64_t,2> res = {x,y};
	return res;
}


std::array<std::vector<std::vector<GRBVar>>,4> 
addSpeckKSRXDConstr(GRBModel & model,
				 uint const n,
				 uint const m,
				 uint const nbRound,
				 uint const alpha,
				 uint const beta,
				 uint const gamma){
//Add constraints for the key schedule of Speck for RX-differentials
//Return the master key, round key, l and tmpKS variables

	//Create variables for the key schedule
	auto mk = genArray2DBinVar(model, m, n, "dmk");
    auto k = genArray2DBinVar(model, nbRound, n, "dk");
	//We need the l[i] variables from the specification, but also some additional temp variables
	auto l = genArray2DBinVar(model, nbRound+m-2, n, "dl");
	auto tmpKS = genArray2DBinVar(model, nbRound-1, n, "dtmpKS");
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
	for(uint64_t i = 0; i < nbRound-1; i++){
		vector<GRBVar> Sali(n); //l[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Sali[j] = l[i][mod(j+alpha,n)];

		vector<GRBVar> Sbki(n); //k[i] rotated to the left by beta
		for(uint j = 0; j < n; j++)
			Sbki[j] = k[i][mod(j-beta,n)];

		//l[i+m-1] = (k[i] + Sali) ^ i
		//First, tmp = k[i] + Sali
		addModAddRXDiffConstr(model, k[i], Sali, tmpKS[i], gamma);

		//l[i+m-1] = tmp ^ i
		//Write i as a binary constant vector
		vector<uint> bini(n);
		for(uint j = 0; j < n; j++)
			bini[j] = (i >> j)&1;
		addRXCstXORConstr(model, tmpKS[i], bini, l[i+m-1], gamma);

		//k[i+1] = Sbki ^ l[i+m-1]
		for(uint j = 0; j < n; j++)
			addXORConstr(model, Sbki[j], l[i+m-1][j], k[i+1][j]);
	}

	array<vector<vector<GRBVar>>,4> ret = {mk,k,l,tmpKS};
    return ret;
}

std::array<std::vector<std::vector<GRBVar>>,3> 
addSpeckRXDConstrRelatedKey(GRBModel & model,
				  uint const n,
				  uint const nbRound,
				  uint const alpha, 
				  uint const beta, 
				  uint const gamma,
				  std::vector<std::vector<GRBVar>> & k){
//Add constraints for RX-differential propagation in Speck in related key
//Returns a triplet for containing the variables [x,y,z]

	//State variables
	auto x = genArray2DBinVar(model, nbRound+1, n, "dx");
	auto y = genArray2DBinVar(model, nbRound+1, n, "dy");
	auto z = genArray2DBinVar(model, nbRound, n, "dz");

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

		for(uint j = 0; j < n; j++){
			//y[i+1] = Sbyi ^ x[i+1]
			addXORConstr(model, Sbyi[j], x[i+1][j], y[i+1][j]);
			//x[i+1] = z[i] ^ k[i]
			addXORConstr(model, z[i][j], k[i][j], x[i+1][j]);
		}
	}

	array<vector<vector<GRBVar>>,3> ret = {x,y,z};
    return ret;
}

GRBModel
getSpeckRelatedKeyModel(uint const n,
			  uint const m, 
			  uint const nbRound,
			  uint const alpha, 
			  uint const beta, 
			  uint const gamma,
			  GRBEnv & env){
//Create a model for related-key RX-diff for Speck 2n/mn over nbRounds rounds, using alpha/beta for the round function. env is only here for better console login from gurobi
//Returns the model

	// Create an empty model
    GRBModel model = GRBModel(env);
    //Key schedule constraints
	auto [mk,k,l,tmpKS] = addSpeckKSRXDConstr(model,n,m,nbRound,alpha,beta,gamma);

	//Constraints for the differential trail
    auto [x,y,z] = addSpeckRXDConstrRelatedKey(model,n,nbRound,alpha,beta,gamma,k);

    //notAllEqual variables
    vector<vector<GRBVar>> neq(nbRound, vector<GRBVar>(n-1));
    vector<vector<GRBVar>> neqk(nbRound-1, vector<GRBVar>(n-1));
    for(uint i = 0; i < nbRound; i++){
    	for(uint j = 0; j < gamma-1; j++)
    		neq[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "neq"+to_string(i)+"_"+to_string(j));
    	//don't check the bit of index gamma-1, as it's the MSB of the right half
    	for(uint j = gamma; j < n-1; j++)
    		neq[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "neq"+to_string(i)+"_"+to_string(j));
    }
    for(uint i = 0; i < nbRound-1; i++){
    	for(uint j = 0; j < gamma-1; j++)
    		neqk[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "neqk"+to_string(i)+"_"+to_string(j));
    	//don't check the bit of index gamma-1, as it's the MSB of the right half
    	for(uint j = gamma; j < n-1; j++)
    		neqk[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "neq"+to_string(i)+"_"+to_string(j));
    }

    //Sum of neq vars
    //Could be useful to know the split between rounds and key schedules, so introduce two intermediary variables for that
    //Shouldn't make the model more complicated as they'll get cleared up by the presolver, but way easier to handle other constraints later on
    
	//sum of all the ones involved in the same state/key 
	vector<GRBVar> sumneq_state(nbRound);
	vector<GRBVar> sumneq_key(nbRound-1); //one less in key schedule
	for(uint r = 0; r < nbRound; r++)
		sumneq_state[r] = model.addVar(0.0, n, 0.0, GRB_INTEGER, "sumneq_state"+to_string(r));
	for(uint r = 0; r < nbRound-1; r++)
		sumneq_key[r] = model.addVar(0.0, n, 0.0, GRB_INTEGER, "sumneq_key"+to_string(r));
	//Also the sum of both of the above involved in a single rounge
	//Note that for round 0, this is just the state (no key scheduling needed)
	//For round i > 0, it's sumneq_state[i] + sumneq_key[i-1]
	vector<GRBVar> sumneq_round(nbRound);
	sumneq_round[0] = model.addVar(0.0, n, 0.0, GRB_INTEGER, "sumneq_round0");
	for(uint r = 1; r < nbRound; r++)
		sumneq_round[r] = model.addVar(0.0, 2*n, 0.0, GRB_INTEGER, "sumneq_round"+to_string(r));
	//Finally, sum over all states, all key, all everything
	GRBVar sumneq_allstate = model.addVar(0.0, nbRound*n, 0.0, GRB_INTEGER, "sumneq_allstate");
    GRBVar sumneq_allkey = model.addVar(0.0, nbRound*n, 0.0, GRB_INTEGER, "sumneq_allkey");
    GRBVar sumneq = model.addVar(0.0, 2*nbRound*n, 0.0, GRB_INTEGER, "sumneq");

    //Sum neq on states
    GRBLinExpr expr_sumneq_allstate = 0;
    for(uint r = 0; r < nbRound; r++){
    	GRBLinExpr expr_sumneq_state = 0;
		for(uint i = 0; i < gamma-1; i++)
			expr_sumneq_state += neq[r][i];
		for(uint i = gamma; i < n-1; i++)
			expr_sumneq_state += neq[r][i];
		model.addConstr(sumneq_state[r] == expr_sumneq_state);
		expr_sumneq_allstate += sumneq_state[r];
	}
	model.addConstr(sumneq_allstate == expr_sumneq_allstate);
	
	//Sum neq on key schedule
	GRBLinExpr expr_sumneq_key = 0;
    for(uint r = 0; r < nbRound-1; r++){
    	GRBLinExpr expr_sumneq_key = 0;
		for(uint i = 0; i < gamma-1; i++)
			expr_sumneq_key += neqk[r][i];
		for(uint i = gamma; i < n-1; i++)
			expr_sumneq_key += neqk[r][i];
		model.addConstr(sumneq_key[r] == expr_sumneq_key);
		expr_sumneq_key += sumneq_key[r];
	}
	model.addConstr(sumneq_allkey == expr_sumneq_key);

	//Sum in the same round
	//Round 0 is just the state
	model.addConstr(sumneq_round[0] == sumneq_state[0]);
	//Otherwise, state r + key r-1
	for(uint r = 1; r < nbRound; r++)
		model.addConstr(sumneq_round[r] == sumneq_state[r] + sumneq_key[r-1]);

	//Finally, sum of everything
	model.addConstr(sumneq == sumneq_allkey + sumneq_allstate);


	//Proper constraints on neq vars to compute them
	//notAllEqual constraints rounds
	for(uint i = 0; i < nbRound; i++){
		vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Saxi[j] = x[i][mod(j+alpha,n)];

		//notAllEq constraints
    	for(uint j = 0; j < gamma-1; j++){
    		vector<GRBVar> vars({Saxi[j], y[i][j], z[i][j]});
    		addNotAllEqualBoolConstr(model, vars, neq[i][j]);
    	}
    	for(uint j = gamma; j < n-1; j++){
    		vector<GRBVar> vars({Saxi[j], y[i][j], z[i][j]});
    		addNotAllEqualBoolConstr(model, vars, neq[i][j]);
    	}
	}

	//notAllEqual constraints key schedule
	for(uint i = 0; i < nbRound-1; i++){
		vector<GRBVar> Sali(n); //l[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Sali[j] = l[i][mod(j+alpha,n)];

		//notAllEq constraints
    	for(uint j = 0; j < gamma-1; j++){
    		vector<GRBVar> vars({k[i][j], Sali[j], tmpKS[i][j]});
    		addNotAllEqualBoolConstr(model, vars, neq[i][j]);
    	}
    	for(uint j = gamma; j < n-1; j++){
    		vector<GRBVar> vars({k[i][j], Sali[j], tmpKS[i][j]});
    		addNotAllEqualBoolConstr(model, vars, neq[i][j]);
    	}
	}
	model.update();
	return model;
}

std::pair<std::vector<GRBVar>, std::vector<GRBVar>>
addModAddProbaSpeckRelatedKey(GRBModel & model,
							  uint const alpha,
							  uint const gamma,
							  std::vector<std::vector<GRBVar>> & x,
							  std::vector<std::vector<GRBVar>> & y,
							  std::vector<std::vector<GRBVar>> & z,
							  std::vector<std::vector<GRBVar>> & k,
							  std::vector<std::vector<GRBVar>> & l,
							  std::vector<std::vector<GRBVar>> & tmpKS){
//Add mod add probabilities constraints and variables for the related key model
//Returns a pair of vector of variables, .first the weight variables for the state, .second for the weight variables of the key schedule
	uint nbRound = x.size()-1;
	uint n = x[0].size();

	//Weight variables
    vector<GRBVar> wtState(nbRound);
    vector<GRBVar> wtKey(nbRound-1);
	for(uint i = 0; i < nbRound; i++)
		wtState[i] = model.addVar(0,2*n,0,GRB_CONTINUOUS, "wtState"+to_string(i));
	for(uint i = 0; i < nbRound-1; i++)
		wtKey[i] = model.addVar(0,2*n,0,GRB_CONTINUOUS, "wtKey"+to_string(i));
	model.update();

	//Mod add in states
	for(uint i = 0; i < nbRound; i++){
		vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Saxi[j] = x[i][mod(j+alpha,n)];

    	//Probability constraints
    	addModAddRXProbaConstr(model,Saxi,y[i],z[i],wtState[i],gamma,"state_r"+to_string(i));
	}

	//Mod add in key schedule
	for(uint i = 0; i < nbRound-1; i++){
		vector<GRBVar> Sali(n); //l[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Sali[j] = l[i][mod(j+alpha,n)];

		addModAddRXProbaConstr(model, k[i], Sali, tmpKS[i], wtKey[i],gamma,"key_r"+to_string(i));
	}
	model.update();

	return make_pair(wtState, wtKey);
}


std::tuple<double,uint, std::vector<std::array<uint64_t,3>>> 
searchBestTrailSpeckRelatedKey(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				bool const withExactProbability,
				std::map<uint,uint> const & lowerBoundConsecutiveRounds,
				std::vector<uint64_t> const & inputDiff,
				std::vector<uint64_t> const & outputDiff){
	uint val_sumneq = 0;
	vector<array<uint64_t,3>> trail;
	auto knownTrail = make_pair(val_sumneq, trail);
	return searchBestTrailSpeckRelatedKey(params,gamma,env,withExactProbability,lowerBoundConsecutiveRounds,knownTrail,inputDiff,outputDiff);
}
std::tuple<double,uint, std::vector<std::array<uint64_t,3>>> 
searchBestTrailSpeckRelatedKey(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				bool const withExactProbability,
				std::map<uint,uint> const & lowerBoundConsecutiveRounds,
				std::pair<uint, std::vector<std::array<uint64_t,3>>> const & knownTrail,
				std::vector<uint64_t> const & inputDiff,
				std::vector<uint64_t> const & outputDiff){
//Same as existTrailSpeck but for related key
//the returned trail/the one used as input will contain the key RXdiff trail instead of key value

	auto [n,m,nbRound,alphaRot,betaRot] = params;
	auto model = getSpeckRelatedKeyModel(n,m,nbRound,alphaRot,betaRot,gamma,env);

	//variables
	auto x = get2DArrayVar(model, nbRound+1, n, "dx");
    auto y = get2DArrayVar(model, nbRound+1, n, "dy");
    auto z = get2DArrayVar(model, nbRound, n, "dz"); 
    auto mk = get2DArrayVar(model, m, n, "dmk");
    auto k = get2DArrayVar(model, nbRound, n, "dk");
    auto l = get2DArrayVar(model, nbRound+m-2, n, "dl");
	auto tmpKS = get2DArrayVar(model, nbRound-1, n, "dtmpKS");
    GRBVar sumneq = model.getVarByName("sumneq");

	//Input vars and constraints (if any)
	if(inputDiff.size() > 0){
		addValueToBinVectorConstr(model, x[0], inputDiff[0]);
		addValueToBinVectorConstr(model, y[0], inputDiff[1]);
	}
	if(outputDiff.size() > 0){
		addValueToBinVectorConstr(model, x[nbRound], outputDiff[0]);
		addValueToBinVectorConstr(model, y[nbRound], outputDiff[1]);
	}

    //Prepare the callback
	vector<array<vector<GRBVar>,3>> modAddVars(nbRound+nbRound-1);
	for(uint r = 0; r < nbRound; r++){
		modAddVars[r][0] = vector<GRBVar>(n);
		modAddVars[r][1] = vector<GRBVar>(n);
		modAddVars[r][2] = vector<GRBVar>(n);
		for(uint i = 0; i < n; i++){
			modAddVars[r][0][i] = x[r][mod(i+alphaRot,n)];
			modAddVars[r][1][i] = y[r][i];
			modAddVars[r][2][i] = z[r][i];
		}
	}
	for(uint r = 0; r < nbRound-1; r++){
		modAddVars[nbRound+r][0] = vector<GRBVar>(n);
		modAddVars[nbRound+r][1] = vector<GRBVar>(n);
		modAddVars[nbRound+r][2] = vector<GRBVar>(n);
		for(uint i = 0; i < n; i++){
			modAddVars[nbRound+r][0][i] = k[r][i];
			modAddVars[nbRound+r][1][i] = l[r][mod(i+alphaRot,n)];
			modAddVars[nbRound+r][2][i] = tmpKS[r][i];
		}
	}
	vector<vector<GRBVar>> stateVars(nbRound+1, vector<GRBVar>(2*n));
	for(uint r = 0; r < nbRound+1; r++){
		for(uint i = 0; i < n; i++){
			stateVars[r][i] = x[r][i];
			stateVars[r][i+n] = y[r][i];
		}
	}
	model.update();
	CustomCallback cb(n,gamma,modAddVars,stateVars,sumneq,k);
	model.setCallback(&cb);

	//Check if we have a valid trail provided as input, if so, use it to bound the number of neq variables and provide a starting value for the trail variables
	auto const & trail = knownTrail.second;
    if(trail.size() > 0){
		//Bound the number of notAllEq
		// auto const bestBound = knownTrail.first;
		// model.addConstr(sumneq <= bestBound);

		//Provide MIP Start info from the trail
		for(uint r = 0; r < nbRound+1; r++){
			auto const trail_r = trail[r];
			for(uint i = 0; i < n; i++){
				uint val = (trail_r[0] >> i)&1;
				x[r][i].set(GRB_DoubleAttr_Start, val);

				val = (trail_r[1] >> i)&1;
				y[r][i].set(GRB_DoubleAttr_Start, val);

				if(r < nbRound){
					val = (trail_r[2] >> i)&1;
					k[r][i].set(GRB_DoubleAttr_Start, val);
				}
			}
		}
    }

    //Lower bounds based on results on a smaller number of rounds
    if(lowerBoundConsecutiveRounds.size() > 0){
    	//Get the variable representing the sum of neq vars over a given round
    	auto sumneq_round = getArrayVar(model,nbRound,"sumneq_round");
    	addMatsuiLikeConstr(model, lowerBoundConsecutiveRounds, sumneq_round, sumneq);
    }

    //Set the objective depending on withExactProbability
    if(withExactProbability){
    	auto [wtState, wtKey] = addModAddProbaSpeckRelatedKey(model,alphaRot,gamma,x,y,z,k,l,tmpKS);
    	GRBLinExpr obj = 0;
		for(uint i = 0; i < nbRound; i++)
			obj += wtState[i];
		for(uint i = 0; i < nbRound-1; i++)
			obj += wtKey[i];
		//Objective: minimize probability
		model.setObjective(obj, GRB_MINIMIZE);
    }
    else{
    	//Objective: minimize the number of notAllEq
		GRBLinExpr obj = sumneq;
		model.setObjective(obj, GRB_MINIMIZE);
    }

    model.update();
	model.optimize();
	if(model.get(GRB_IntAttr_SolCount) > 0){

	    vector<uint64_t> valx(nbRound+1,0);
	    vector<uint64_t> valy(nbRound+1,0);
	    vector<uint64_t> valz(nbRound,0);
	    vector<uint64_t> valk(nbRound,0);
	    vector<uint64_t> vall(nbRound+m-2);
	    vector<uint64_t> valtmpKS(nbRound-1);
	    vector<uint64_t> valmk(m,0);
	    vector<uint64_t> alphaState(nbRound,0);
	    vector<uint64_t> betaState(nbRound,0);
	    vector<uint64_t> deltaState(nbRound,0);
	    vector<uint64_t> alphaKey(nbRound-1,0);
	    vector<uint64_t> betaKey(nbRound-1,0);
	    vector<uint64_t> deltaKey(nbRound-1,0);
	    vector<double> modAddLogState(nbRound,0);
	    vector<double> modAddLogKey(nbRound-1,0);
	    double sumLogState = 0;
	    double sumLogKey = 0;

	    //Get the values of variables in the solution
	    for(uint i = 0; i < nbRound+1; i++){
	    	valx[i] = getUintSolutionFromBinVector(x[i]);
	    	valy[i] = getUintSolutionFromBinVector(y[i]);
	    }
	    for(uint i = 0; i < nbRound; i++){
	    	valz[i] = getUintSolutionFromBinVector(z[i]);
	    	valk[i] = getUintSolutionFromBinVector(k[i]);
	    }
	    for(uint i = 0; i < nbRound+m-2; i++)
	    	vall[i] = getUintSolutionFromBinVector(l[i]);
	    for(uint i = 0; i < nbRound-1; i++)
	    	valtmpKS[i] = getUintSolutionFromBinVector(tmpKS[i]);

	    for(uint i = 0; i < m; i++)
			valmk[i] = getUintSolutionFromBinVector(mk[i]);

		//Get the mod add values in the state
	    for(uint i = 0; i < nbRound; i++){
	    	//zi = (xi >>> alphaRot) + yi
	    	vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
			for(uint j = 0; j < n; j++)
				Saxi[j] = x[i][mod(j+alphaRot,n)];
			uint64_t tmpint = getUintSolutionFromBinVector(Saxi);

			//Keep the differences used in modular additions
	    	alphaState[i] = tmpint;
	    	betaState[i] = valy[i];
	    	deltaState[i] = valz[i];

	    	uint64_t countsol = getRXDiffCount(alphaState[i],betaState[i],deltaState[i],n,gamma);
			sumLogState += (log2(countsol)-2*n);
			modAddLogState[i] = log2(countsol)-2*n;
		}

		//Get the mod add values in the key schedule
		for(uint i = 0; i < nbRound-1; i++){
			vector<GRBVar> Sali(n); //l[i] rotated to the right by alpha
			for(uint j = 0; j < n; j++)
				Sali[j] = l[i][mod(j+alphaRot,n)];
			uint64_t tmpint = getUintSolutionFromBinVector(Sali);

			//Keep the difference used in modular additions
			alphaKey[i] = valk[i];
			betaKey[i] = tmpint;
			deltaKey[i] = valtmpKS[i];

			uint64_t countsol = getRXDiffCount(alphaKey[i], betaKey[i], deltaKey[i],n,gamma);
			sumLogKey += (log2(countsol)-2*n);
			modAddLogKey[i] = log2(countsol)-2*n;
		}

		

		uint val_sumneq = uint(round(sumneq.get(GRB_DoubleAttr_X)));
		auto bestObj = model.get(GRB_DoubleAttr_ObjVal);
		auto bestLB = model.get(GRB_DoubleAttr_ObjBound);

		cout << "***********************************************************************" << endl;
		cout << "*** Speck " << 2*n << "/" << m*n << " " << nbRound << " rounds with k = " << gamma << " related key ***" << endl;
		cout << "*** Best solution found with objective " << bestObj << " ***" << endl;
		if(model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
			cout << "Solution obtained after time-out" << endl;
		else
			cout << "Solution obtained after full optimization" << endl;

		cout << "Number of notAllEqual vars : " << val_sumneq << endl;
		cout << "Best lower bound on objective : " << bestLB << endl;

		cout << "Trail :" << endl;
		cout << "Master key : ";
		for(uint i = 0; i < m; i++)
			cout << "0x" << setfill('0') << setw(n/4) << hex << valmk[i] << " ";
		cout << endl;
		for(uint i = 0; i < nbRound+1; i++){
			cout <<  "x" << i << " : " << "0x" << setfill('0') << setw(n/4) << hex << valx[i];
			cout << " y" << i << " : " << "0x" << setfill('0') << setw(n/4) << hex << valy[i];
			if(i < nbRound)
				cout << " k" << i << " : " << "0x" << setfill('0') << setw(n/4) << hex << valk[i];
			cout << endl;
		}
		cout << "modAdd State:" << endl;
		for(uint i = 0; i < nbRound; i++){
			cout << "0x" << setfill('0') << setw(n/4) << hex << alphaState[i] << " + ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << betaState[i] << " -> ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << deltaState[i] << dec;
			cout << " (2^" << modAddLogState[i] << ")" << endl;
		}
		cout << "sumLogState = " << sumLogState << endl;
		cout << "modAddTrail = [";
		for(uint i = 0; i < nbRound; i++){
			cout << "[";
			cout << "0x" << setfill('0') << setw(n/4) << hex << alphaState[i];
			cout << ",";
			cout << "0x" << setfill('0') << setw(n/4) << hex << betaState[i];
			cout << ",";
			cout << "0x" << setfill('0') << setw(n/4) << hex << deltaState[i];
			if(i < nbRound-1)
				cout << "]," << endl;
			else
				cout << "]]" << endl;
		}
		cout << "modAdd Key:" << endl;
		for(uint i = 0; i < nbRound-1; i++){
			cout << "0x" << setfill('0') << setw(n/4) << hex << alphaKey[i] << " + ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << betaKey[i] << " -> ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << deltaKey[i] << dec;
			cout << " (2^" << modAddLogKey[i] << ")" << endl;
		}
		cout << "sumLogKey = " << sumLogKey << endl;
		cout << "modAddTrail = [";
		for(uint i = 0; i < nbRound-1; i++){
			cout << "[";
			cout << "0x" << setfill('0') << setw(n/4) << hex << alphaKey[i];
			cout << ",";
			cout << "0x" << setfill('0') << setw(n/4) << hex << betaKey[i];
			cout << ",";
			cout << "0x" << setfill('0') << setw(n/4) << hex << deltaKey[i];
			if(i < nbRound-1)
				cout << "]," << endl;
			else
				cout << "]]" << endl;
		}
		cout << "***********************************************************************" << endl;
		cout << dec;


		vector<array<uint64_t,3>> trail(nbRound+1);
		for(uint r = 0; r < nbRound+1; r++){
			trail[r][0] = valx[r];
			trail[r][1] = valy[r];
			if(r < nbRound)
				trail[r][2] = valk[r];
		}

		return make_tuple(bestLB, val_sumneq, trail);
	}
	else{
		cout << "No solution" << endl;
		uint val_sumneq = 0;
		double bestLB = 0;
		vector<array<uint64_t,3>> trail;
		return make_tuple(bestLB, val_sumneq, trail);
	}
}

GRBModel
getSpeckModelPlaintextPairRelatedKey(uint const n,
						   uint const m, 
						   uint const nbRound,
						   uint const alpha, 
						   uint const beta, 
						   uint const gamma,
						   GRBEnv & env,
						   std::vector<std::array<uint64_t,3>> const & trail){
//Create a MILP model for Speck to find a plaintext pair matching a given trail in related key

	// Create an empty model
    GRBModel model = GRBModel(env);

    //Key schedule constraints
	auto [mk0,k0] = addSpeckKSConstr(model,n,m,nbRound,alpha,beta,"0");
	auto [mk1,k1] = addSpeckKSConstr(model,n,m,nbRound,alpha,beta,"1");

	//Constraints for the differential trail
	auto [mk,k,l,tmpKS] = addSpeckKSRXDConstr(model,n,m,nbRound,alpha,beta,gamma);
    auto [x,y,z] = addSpeckRXDConstrRelatedKey(model,n,nbRound,alpha,beta,gamma,k);

    //Constraints for the plaintext values
    auto [x0,y0,z0] = addSpeckPlaintextValueConstr(model,n,nbRound,alpha,beta,"0",k0);
    auto [x1,y1,z1] = addSpeckPlaintextValueConstr(model,n,nbRound,alpha,beta,"1",k1);

    //--- Bind the value variables to the differential variables ---
    //Essentially, enforce dx = rot(x0) xor x1
    for(uint r = 0; r < nbRound+1; r++){
    	addRXDifferentialValueBinding(model,x0[r],x1[r],x[r],gamma);
    	addRXDifferentialValueBinding(model,y0[r],y1[r],y[r],gamma);
    }
    for(uint r = 0; r < nbRound; r++){
    	addRXDifferentialValueBinding(model,z0[r],z1[r],z[r],gamma);
    	addRXDifferentialValueBinding(model,k0[r],k1[r],k[r],gamma);
    }
    for(uint r = 0; r < m; r++)
    	addRXDifferentialValueBinding(model,mk0[r],mk1[r],mk[r],gamma);

    //Add the constraint values for the trail
    //Fixing the value of (x,y,k) is enough to deduce a unique value for z
    //The solver should handle that easily during presolving
    for(uint r = 0; r < nbRound+1; r++){
    	addValueToBinVectorConstr(model,x[r],trail[r][0]);
    	addValueToBinVectorConstr(model,y[r],trail[r][1]);
    }
    for(uint r = 0; r < nbRound; r++)
    	addValueToBinVectorConstr(model,k[r],trail[r][2]);

    //Some parameters that should help solving
    //We're only interested in getting one solution
    model.set(GRB_IntParam_SolutionLimit, 1);
    //Focus on finding solutions rather than optimization
    model.set(GRB_IntParam_MIPFocus, 1);
    //Still give an arbitrary objective as it's supposed to help the solver a bit
    GRBLinExpr obj = 0;
    for(uint i = 0; i < n; i++){
    	obj += x0[0][i];
    	obj += y0[0][i];
    }
    model.setObjective(obj, GRB_MINIMIZE);

    return model;
}

bool
findPlaintextForTrailSpeckRelatedKey(std::array<uint,5> const & params,
				uint const gamma,
				GRBEnv & env,
				std::vector<std::array<uint64_t,3>> const & trail){

	auto [n,m,nbRound,alpha,beta] = params;
	auto model = getSpeckModelPlaintextPairRelatedKey(n,m,nbRound,alpha,beta,gamma,env,trail);
	model.optimize();

	//Extract the solution, if any
	if(model.get(GRB_IntAttr_SolCount) == 0){
		cout << "*************************************" << endl;
		cout << "*** No matching pair of plaintext ***" << endl;
		cout << "*************************************" << endl;
		return false;
	}
	else{

		auto x0 = getArrayVar(model,n,"x00_");
		auto x1 = getArrayVar(model,n,"x10_");
		auto y0 = getArrayVar(model,n,"y00_");
		auto y1 = getArrayVar(model,n,"y10_");
		auto k = get2DArrayVar(model,nbRound,n,"dk");
		auto k0 = get2DArrayVar(model,nbRound,n,"k0");
		auto k1 = get2DArrayVar(model,nbRound,n,"k1");
		auto mk = get2DArrayVar(model,m,n,"dmk");
		auto mk0 = get2DArrayVar(model,m,n,"mk0");
		auto mk1 = get2DArrayVar(model,m,n,"mk1");

		array<uint64_t,4> res;
		res[0] = getUintSolutionFromBinVector(x0);
		res[1] = getUintSolutionFromBinVector(y0);
		res[2] = getUintSolutionFromBinVector(x1);
		res[3] = getUintSolutionFromBinVector(y1);
		vector<uint64_t> valk(nbRound);
		vector<uint64_t> valk0(nbRound);
		vector<uint64_t> valk1(nbRound);
		vector<uint64_t> valmk(m);
		vector<uint64_t> valmk0(m);
		vector<uint64_t> valmk1(m);
		vector<uint64_t> dmk(m);

		for(uint r = 0; r < nbRound; r++){
			valk[r] = getUintSolutionFromBinVector(k[r]);
			valk0[r] = getUintSolutionFromBinVector(k0[r]);
			valk1[r] = getUintSolutionFromBinVector(k1[r]);
		}
		for(uint i = 0; i < m; i++){
			valmk[i] = getUintSolutionFromBinVector(mk[i]);
			valmk0[i] = getUintSolutionFromBinVector(mk0[i]);
			valmk1[i] = getUintSolutionFromBinVector(mk1[i]);
			dmk[i] = CSHL(valmk0[i],gamma,n) ^ valmk1[i];
		}

		uint64_t dx = CSHL(res[0],gamma,n) ^ res[2];
		uint64_t dy = CSHL(res[1],gamma,n) ^ res[3];
		auto [cx0,cy0] = speckEncrypt(n,nbRound,alpha,beta,res[0],res[1],valk0);
		auto [cx1,cy1] = speckEncrypt(n,nbRound,alpha,beta,res[2],res[3],valk1);
		uint64_t cdx = CSHL(cx0,gamma,n) ^ cx1;
		uint64_t cdy = CSHL(cy0,gamma,n) ^ cy1;

		cout << "*******************************" << endl;
		cout << "*** Found pair of plaintext ***" << endl;
		cout << "x0 = 0x" << setfill('0') << setw(n/4) << hex << res[0] << " ";
		cout << "y0 = 0x" << setfill('0') << setw(n/4) << hex << res[1] << " ";
		cout << "mk0 = ";
		for(auto const tmp : valmk0)
			cout << "0x" << setfill('0') << setw(n/4) << hex << tmp << " ";
		cout << endl;
		cout << "x1 = 0x" << setfill('0') << setw(n/4) << hex << res[2] << " ";
		cout << "y1 = 0x" << setfill('0') << setw(n/4) << hex << res[3] << " ";
		cout << "mk1 = ";
		for(auto const tmp : valmk1)
			cout << "0x" << setfill('0') << setw(n/4) << hex << tmp << " ";
		cout << endl;

		//Check the found solution using the actual speck encryption
		cout << "Resulting input dx dy : ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << dx << " ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << dy << endl;
		cout << " Expected input dx dy : ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << trail[0][0] << " ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << trail[0][1] << endl;

		cout << "Resulting output dx dy : ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << cdx << " ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << cdy << endl;
		cout << " Expected output dx dy : ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << trail[nbRound][0] << " ";
		cout << "0x" << setfill('0') << setw(n/4) << hex << trail[nbRound][1] << endl;

		cout << "Resulting master key difference : ";
		for(auto const tmp : dmk)
			cout << "0x" << setfill('0') << setw(n/4) << hex << tmp << " ";
		cout << endl;
		cout << "Expected master key difference : ";
		for(auto const tmp : valmk)
			cout << "0x" << setfill('0') << setw(n/4) << hex << tmp << " ";
		cout << endl;

		auto expk0 = speckKS(n,nbRound,alpha,beta,valmk0);
		auto expk1 = speckKS(n,nbRound,alpha,beta,valmk1);
		cout << "Key schedule : " << endl;
		for(uint r = 0; r < nbRound; r++){
			cout << "k0[" << r << "] = " << "0x" << setfill('0') << setw(n/4) << hex << valk0[r] << " ";
			cout << "k1[" << r << "] = " << "0x" << setfill('0') << setw(n/4) << hex << valk1[r] << " ";
			cout << "expk0[" << r << "] = " << "0x" << setfill('0') << setw(n/4) << hex << expk0[r] << " ";
			cout << "expk1[" << r << "] = " << "0x" << setfill('0') << setw(n/4) << hex << expk1[r] << " ";
			cout << endl;
		}

		if(dx == trail[0][0] && dy == trail[0][1])
			cout << "Input differentials matches" << endl;
		else
			cout << "Input differentials fail" << endl;
		if(cdx == trail[nbRound][0] && cdy == trail[nbRound][1])
			cout << "Output differentials matches" << endl;
		else
			cout << "Output differentials fail" << endl;

		cout << "*******************************" << endl;

		cout << dec;
		return true;
	}
}
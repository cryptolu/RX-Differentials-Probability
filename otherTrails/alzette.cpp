#include "alzette.hpp"
using namespace std;


std::array<std::vector<std::vector<GRBVar>>,3> 
addAlzetteRXDConstr(GRBModel & model,
					std::vector<uint> const & params,
					uint64_t const cst,
					uint const gamma){
//Add constraints for RX-differential propagation in Alzette
//Returns a triplet for containing the variables [x,y,z]
	uint nbRound = params.size()/2;
	uint n = 32;
	vector<uint> bincst(n);
	for(uint i = 0; i < n; i++)
		bincst[i] = (cst >> i)&1;

	//State variables
	auto x = genArray2DBinVar(model, nbRound+1, n, "dx");
	auto y = genArray2DBinVar(model, nbRound+1, n, "dy");
	auto z = genArray2DBinVar(model, nbRound, n, "dz");

    for(uint r = 0;  r < nbRound; r++){
    	vector<GRBVar> tmp(n);
    	//z[r] = x[r] + (y[r] >>> params[2r])
    	for(uint i = 0; i < n; i++)
    		tmp[i] = y[r][(i+params[2*r])%n]; //tmp = y[r] >>> params[2r]
    	//Validity constraint
    	addModAddRXDiffConstr(model, x[r], tmp, z[r], gamma);

    	//x[r+1] = z[r] ^ cst
    	addRXCstXORConstr(model, z[r], bincst, x[r+1], gamma);

    	//y[r+1] = y[r] ^ (z[r] >>> params[2r+1])
    	for(uint i = 0; i < n; i++)
    		tmp[i] = z[r][(i+params[2*r+1])%n]; //tmp = z[r] >>> params[2r+1]
    	for(uint i = 0; i < n; i++)
    		addXORConstr(model, y[r][i], tmp[i], y[r+1][i]);
    }
    array<vector<vector<GRBVar>>,3> ret = {x,y,z};
    return ret;
}

GRBModel getAlzetteModel(std::vector<uint> const & params,
						 uint64_t const cst,
						 uint const gamma,
						 GRBEnv & env){

	uint nbRound = params.size()/2;
	uint n = 32;

	// Create an empty model
    GRBModel model = GRBModel(env);

    //Constraints for the differential trail
    auto [x,y,z] = addAlzetteRXDConstr(model,params,cst,gamma);

    //notAllEqual variables
    vector<vector<GRBVar>> v(nbRound, vector<GRBVar>(n-1));
    for(uint i = 0; i < nbRound; i++){
    	for(uint j = 0; j < gamma-1; j++)
    		v[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "neq"+to_string(i)+"_"+to_string(j));
    	//don't check the bit of index gamma-1, as it's the MSB of the right half
    	for(uint j = gamma; j < n-1; j++)
    		v[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "neq"+to_string(i)+"_"+to_string(j));
    }
    //Sum of all neq vars
    GRBVar sumneq = model.addVar(0.0, nbRound*n, 0.0, GRB_INTEGER, "sumneq");
    model.update();

    GRBLinExpr expr_sumneq = 0;
	for(uint r = 0; r < nbRound; r++){
		for(uint i = 0; i < gamma-1; i++)
			expr_sumneq += v[r][i];
		for(uint i = gamma; i < 31; i++)
			expr_sumneq += v[r][i];
	}
	model.addConstr(expr_sumneq == sumneq);

    model.update();

    for(uint r = 0;  r < nbRound; r++){
    	vector<GRBVar> tmp(32);
    	//z[r] = x[r] + (y[r] >>> params[2r])
    	for(uint i = 0; i < 32; i++)
    		tmp[i] = y[r][(i+params[2*r])%32]; //tmp = y[r] >>> params[2r]

    	//notAllEq constraints
    	for(uint i = 0; i < gamma-1; i++){
    		vector<GRBVar> vars({x[r][i], tmp[i], z[r][i]});
    		addNotAllEqualBoolConstr(model, vars, v[r][i]);
    	}
    	for(uint i = gamma; i < 31; i++){
    		vector<GRBVar> vars({x[r][i], tmp[i], z[r][i]});
    		addNotAllEqualBoolConstr(model, vars, v[r][i]);
    	}
    }
    model.update();

    return model;
}

std::vector<GRBVar> 
addModAddProbaAlzette(GRBModel & model,
					  std::vector<uint> const & params,
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

    for(uint r = 0;  r < nbRound; r++){
    	vector<GRBVar> tmp(n);
    	//z[r] = x[r] + (y[r] >>> params[2r])
    	for(uint i = 0; i < n; i++)
    		tmp[i] = y[r][(i+params[2*r])%n]; //tmp = y[r] >>> params[2r]
    	//Probablity constraint
    	addModAddRXProbaConstr(model, x[r], tmp, z[r], wt[r], gamma, "_r"+to_string(r));
    }
    model.update();

    return wt;
}

std::pair<uint, std::vector<std::array<uint64_t,2>>> 
searchBestTrailAlzette(std::vector<uint> const & params,
				  uint64_t const cst,
				  uint const gamma,
				  GRBEnv & env,
				  bool const withExactProbability,
				  std::vector<uint64_t> const & inputDiff,
				  std::vector<uint64_t> const & outputDiff){

	uint val_sumneq = 0;
	vector<array<uint64_t,2>> trail;
	auto knownTrail = make_pair(val_sumneq, trail);
	return searchBestTrailAlzette(params,cst,gamma,env,withExactProbability,knownTrail,inputDiff,outputDiff);
}

std::pair<int, std::vector<std::array<uint64_t,2>>> 
searchBestTrailAlzette(std::vector<uint> const & params,
				  uint64_t const cst,
				  uint const gamma,
				  GRBEnv & env,
				  bool const withExactProbability,
				  std::pair<uint, std::vector<std::array<uint64_t,2>>> const & knownTrail,
				  std::vector<uint64_t> const & inputDiff,
				  std::vector<uint64_t> const & outputDiff){
/*
Search for the best RX-diff trail with RX rotation gamma
params should be the rotation values  (default {31,24,17,17,0,31,24,16}), number of rounds is deduced for this as params.size()/2
withExactProbability allows to control if we add probability constraints and optimization. False means we just optimize on the number of neq variables
knownTrail can be provided as a pair (W,trail) where W is a bound on the number of neq variables and trail is a starting solution for the model
inputDiff and outputDiff can be left empty to allow any input/output differential, or given to fix the input/output differential respectively, as {x,y} for the difference on the left,right branch. Won't work for n > 64 since differences are given as a uint64_t for each word (but shouldn't really be a problem).
*/

	
	uint n = 32;
	uint nbRound = params.size()/2;
	auto model = getAlzetteModel(params,cst,gamma,env);

	//State variables
    auto x = get2DArrayVar(model, nbRound+1, n, "dx");
    auto y = get2DArrayVar(model, nbRound+1, n, "dy");
    auto z = get2DArrayVar(model, nbRound, n, "dz");

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
			modAddVars[r][0][i] = x[r][i];
			modAddVars[r][1][i] = y[r][(i+params[2*r])%n];
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
	CustomCallback cb(n,gamma,modAddVars,stateVars,sumneq);
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
			}
		}
    }

    //Set the objective depending on withExactProbability
    if(withExactProbability){
    	auto wt = addModAddProbaAlzette(model,params,gamma,x,y,z);
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

			//z[r] = x[r] + (y[r] >>> params[2r])
			//(y[r] >>> params[2r])
	    	vector<GRBVar> tmp(n);
	    	for(uint j = 0; j < n; j++)
	    		tmp[j] = y[i][(j+params[2*i])%n]; //tmp = y[r] >>> params[2r]
	    	uint64_t tmpint = getUintSolutionFromBinVector(tmp);

	    	alpha[i] = valx[i];
	    	beta[i] = tmpint;
	    	delta[i] = valz[i];

	    	uint64_t countsol = getRXDiffCount(alpha[i],beta[i],delta[i],n,gamma);
			sumLog += (log2(countsol)-2*n);
			modAddLog[i] = log2(countsol)-2*n;
		}

		valx[nbRound] = getUintSolutionFromBinVector(x[nbRound]);
	    valy[nbRound] = getUintSolutionFromBinVector(y[nbRound]);

		uint val_sumneq = uint(round(sumneq.get(GRB_DoubleAttr_X)));
		auto bestObj = model.get(GRB_DoubleAttr_ObjVal);
		cout << "***********************************************************************" << endl;
		cout << "*** Alzette 4 rounds with k = " << gamma;
		cout << " cst = " << "0x" << setfill('0') << setw(n/4) << hex << cst << dec;
		cout << " ***" << endl;
		cout << "*** Best solution found with objective " << bestObj << " ***" << endl;
		if(model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
			cout << "Solution obtained after time-out" << endl;
		else
			cout << "Solution obtained after full optimization" << endl;

		cout << "Number of notAllEqual vars : " << val_sumneq << endl;

		cout << "Trail :" << endl;
		for(uint i = 0; i < nbRound+1; i++){
			cout <<  "x" << i << " : " << "0x" << setfill('0') << setw(n/4) << hex << valx[i];
			cout << " y" << i << " : " << "0x" << setfill('0') << setw(n/4) << hex << valy[i];
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
		//std::pair<uint, std::vector<std::array<uint64_t,2>>>

		vector<array<uint64_t,2>> trail(nbRound+1);
		for(uint r = 0; r < nbRound+1; r++){
			trail[r][0] = valx[r];
			trail[r][1] = valy[r];
		}

		return make_pair(val_sumneq, trail);
	}
	else{
		cout << "No solution" << endl;
		uint val_sumneq = 0;
		vector<array<uint64_t,2>> trail;
		return make_pair(val_sumneq, trail);
	}
}


std::pair<uint32_t,uint32_t> 
alzetteCore(std::vector<uint> const & params,
			uint32_t const cst,
			uint32_t const inx,
			uint32_t const iny){
//Alzette core evaluation
	uint nbRound = params.size()/2;
	uint32_t x = inx;
	uint32_t y = iny;
	for(uint r = 0; r < nbRound; r++){
		x += CSHR(y,params[2*r],32);
		y ^= CSHR(x,params[2*r+1],32);
		x ^= cst;	
	}
	return make_pair(x,y);
}

std::array<std::vector<std::vector<GRBVar>>,3> 
addAlzettePlaintextValueConstr(GRBModel & model,
							   std::vector<uint> const & params,
							   uint64_t const cst,
							   std::string const & suffix){
//Add constraints for one plaintext in value for Alzette
// suffix is used to create unique variable names
//Returns a triplet for containing the variables [x,y,z]

	uint n = 32;
	uint nbRound = params.size()/2;
	//Write the constant in binary
	vector<uint> bincst(n);
	for(uint i = 0; i < n; i++)
		bincst[i] = (cst >> i)&1;

	auto x = genArray2DBinVar(model, nbRound+1, n, "x"+suffix);
	auto y = genArray2DBinVar(model, nbRound+1, n, "y"+suffix);
	auto z = genArray2DBinVar(model, nbRound, n, "z"+suffix);
    model.update();

    //Round function constraints
    for(uint r = 0;  r < nbRound; r++){
    	vector<GRBVar> tmp(n);
    	//z[r] = x[r] + (y[r] >>> params[2r])
    	for(uint i = 0; i < n; i++)
    		tmp[i] = y[r][(i+params[2*r])%n]; //tmp = y[r] >>> params[2r]
    	addModAddValueConstr(model, x[r], tmp, z[r], "carry"+suffix+to_string(r));

    	//x[r+1] = z[r] ^ cst
    	for(uint i = 0; i < n; i++){
    		if(bincst[i] == 0)
    			model.addConstr(x[r+1][i] == z[r][i]);
    		else
    			model.addConstr(x[r+1][i] == (1 - z[r][i]));
    	}

    	//y[r+1] = y[r] ^ (z[r] >>> params[2r+1])
    	for(uint i = 0; i < n; i++)
    		tmp[i] = z[r][(i+params[2*r+1])%n]; //tmp = z[r] >>> params[2r+1]
    	for(uint i = 0; i < n; i++)
    		addXORConstr(model, y[r][i], tmp[i], y[r+1][i]);
    }

    array<vector<vector<GRBVar>>,3> ret = {x,y,z};
    return ret;
}

GRBModel
getAlzetteModelPlaintextPair(std::vector<uint> const & params,
							 uint64_t const cst,
							 uint const gamma,
							 GRBEnv & env,
							 std::vector<std::array<uint64_t,2>> const & trail){
//Create a MILP model for Alzette to find a plaintext pair matching a given trail

	uint nbRound = params.size()/2;
	uint n = 32;
	//Write the constant in binary
	vector<uint> bincst(n);
	for(uint i = 0; i < n; i++)
		bincst[i] = (cst >> i)&1;

	// Create an empty model
    GRBModel model = GRBModel(env);

    //--- Differential variable and constraints ---
    auto [x,y,z] = addAlzetteRXDConstr(model,params,cst,gamma);

    //--- Value variables and constraints ---
    auto [x0,y0,z0] = addAlzettePlaintextValueConstr(model,params,cst,"0");
    auto [x1,y1,z1] = addAlzettePlaintextValueConstr(model,params,cst,"1");
    model.update();
    
    //--- Bind the value variables to the differential variables ---
    //Essentially, enforce dx = rot(x0) xor x1
    for(uint r = 0; r < nbRound+1; r++){
    	addRXDifferentialValueBinding(model,x0[r],x1[r],x[r],gamma);
    	addRXDifferentialValueBinding(model,y0[r],y1[r],y[r],gamma);
    }
    for(uint r = 0; r < nbRound; r++)
    	addRXDifferentialValueBinding(model,z0[r],z1[r],z[r],gamma);

    //Add the constraint values for the trail
    //Fixing the value of (x,y) is enough to deduce a unique value for z
    //The solver should handle that easily during presolving
    for(uint r = 0; r < nbRound+1; r++){
    	addValueToBinVectorConstr(model,x[r],trail[r][0]);
    	addValueToBinVectorConstr(model,y[r],trail[r][1]);
    }

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
findPlaintextForTrailAlzette(std::vector<uint> const & params,
				  uint64_t const cst,
				  uint const gamma,
				  GRBEnv & env,
				  std::vector<std::array<uint64_t,2>> const & trail){
//Return a pair {b, p}, b = true if a plaintext pair was found (and then p contains the pair as x0,y0,x1,y1)

	uint n = 32;
	uint nbRound = params.size()/2;
	auto model = getAlzetteModelPlaintextPair(params,cst,gamma,env,trail);
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

		uint64_t dx = CSHL(res[0],gamma,n) ^ res[2];
		uint64_t dy = CSHL(res[1],gamma,n) ^ res[3];
		auto [cx0,cy0] = alzetteCore(params,cst,res[0],res[1]);
		auto [cx1,cy1] = alzetteCore(params,cst,res[2],res[3]);
		uint64_t cdx = CSHL(cx0,gamma,n) ^ cx1;
		uint64_t cdy = CSHL(cy0,gamma,n) ^ cy1;

		cout << "*******************************" << endl;
		cout << "*** Found pair of plaintext ***" << endl;
		cout << "x0 = 0x" << setfill('0') << setw(n/4) << hex << res[0] << " ";
		cout << "y0 = 0x" << setfill('0') << setw(n/4) << hex << res[1] << endl;
		cout << "x1 = 0x" << setfill('0') << setw(n/4) << hex << res[2] << " ";
		cout << "y1 = 0x" << setfill('0') << setw(n/4) << hex << res[3] << endl;

		//Check the found solution using alzetteCore
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

		return make_pair(true,res);
	}

}
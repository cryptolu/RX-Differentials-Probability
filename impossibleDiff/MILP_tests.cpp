#include "MILP_tests.hpp"
using namespace std;

bool testXOR2(bool const dummyVar, uint const c){
	if(c == 0)
		cout << "Test x XOR y = z";
	else
		cout << "Test x XOR y XOR 1 = z";
	if(dummyVar)
		cout << " with dummy variable" << endl;
	else
		cout << " without dummy variable" << endl;
	try {

    // Create an environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
    GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
    GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

    if(dummyVar)
    	addXORConstr(model,x,y,z,c,"dum");	
    else
    	addXORConstr(model,x,y,z,c);

    model.set(GRB_IntParam_PoolSolutions, 1024);
    model.set(GRB_IntParam_PoolSearchMode, 2);
    model.optimize();

    // Print number of solutions stored
    int nSolutions = model.get(GRB_IntAttr_SolCount);
    cout << "Number of solutions found: " << nSolutions;
    if(nSolutions == 4)
    	cout << " : correct number of solutions" << endl;
    else{
    	cout << " : incorrect number of solutions, expected 4" << endl;
    	return false;
    }

    // Print objective values of solutions
    for (int e = 0; e < nSolutions; e++) {
      model.set(GRB_IntParam_SolutionNumber, e);
      uint xv = uint(round(x.get(GRB_DoubleAttr_Xn)));
      uint yv = uint(round(y.get(GRB_DoubleAttr_Xn)));
      uint zv = uint(round(z.get(GRB_DoubleAttr_Xn)));
      cout << "Sol ";
      if(e < 10) cout << " ";
      cout << e << " : ";
      cout << "x=" << xv << " ";
      cout << "y=" << yv << " ";
      cout << "z=" << zv << " ";
      if((xv^yv^c) == zv)
      	cout << "Correct value";
      else{
      	cout << "Incorrect value";
      	return false;
      }
      cout << endl;
    }
    cout << endl;


  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return true;
}

bool testXOR3(bool const dummyVar, uint const c){
	if(c == 0)
		cout << "Test x XOR y XOR w = z";
	else
		cout << "Test x XOR y XOR w XOR 1 = z";
	if(dummyVar)
		cout << " with dummy variable" << endl;
	else
		cout << " without dummy variable" << endl;
	try {

    // Create an environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
    GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
    GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");
    GRBVar w = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "w");

    if(dummyVar)
    	addXORConstr(model,x,y,w,z,c,"dum");	
    else
    	addXORConstr(model,x,y,w,z,c);

    model.set(GRB_IntParam_PoolSolutions, 1024);
    model.set(GRB_IntParam_PoolSearchMode, 2);
    model.optimize();

    // Print number of solutions stored
    int nSolutions = model.get(GRB_IntAttr_SolCount);
    cout << "Number of solutions found: " << nSolutions;
    if(nSolutions == 8)
    	cout << " : correct number of solutions" << endl;
    else{
    	cout << " : incorrect number of solutions, expected 8" << endl;
    	return false;
    }

    // Print objective values of solutions
    for (int e = 0; e < nSolutions; e++) {
      model.set(GRB_IntParam_SolutionNumber, e);
      uint xv = uint(round(x.get(GRB_DoubleAttr_Xn)));
      uint yv = uint(round(y.get(GRB_DoubleAttr_Xn)));
      uint zv = uint(round(z.get(GRB_DoubleAttr_Xn)));
      uint wv = uint(round(w.get(GRB_DoubleAttr_Xn)));
      cout << "Sol ";
      if(e < 10) cout << " ";
      cout << e << " : ";
      cout << "x=" << xv << " ";
      cout << "y=" << yv << " ";
      cout << "w=" << wv << " ";
      cout << "z=" << zv << " ";
      if((xv^yv^wv^c) == zv)
      	cout << "Correct value";
      else{
      	cout << "Incorrect value";
      	return false;
      }
      cout << endl;
    }
    cout << endl;


  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return true;
}

bool testXORn(uint const nvar, uint const c){
	if(c == 0)
		cout << "Test x0 XOR ... XOR x" << nvar-1 << " = z" << endl;
	else
		cout << "Test x0 XOR ... XOR x" << nvar-1 << " + 1 = z" << endl;

	try {

    // Create an environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    vector<GRBVar> vars(nvar);
    for(uint i = 0; i < nvar; i++)
    	vars[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i));
    GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

    addXORConstr(model,vars,z,c,"dum");

    model.set(GRB_IntParam_PoolSolutions, (1ULL << (nvar+1)));
    model.set(GRB_IntParam_PoolSearchMode, 2);
    model.optimize();

    // Print number of solutions stored
    uint nSolutions = model.get(GRB_IntAttr_SolCount);
    cout << "Number of solutions found: " << nSolutions;
    if(nSolutions == (1ULL << nvar))
    	cout << " : correct number of solutions" << endl;
    else{
    	cout << " : incorrect number of solutions, expected " << (1ULL << nvar) << endl;
    	return false;
    }

    for (uint e = 0; e < nSolutions; e++) {
      model.set(GRB_IntParam_SolutionNumber, e);

      vector<uint> vals(nvar);
      for(uint i = 0; i < nvar; i++)
      	vals[i] = uint(round(vars[i].get(GRB_DoubleAttr_Xn)));

      uint zv = uint(round(z.get(GRB_DoubleAttr_Xn)));
    
      cout << "Sol ";
      if(e < 10) cout << " ";
      cout << e << " : ";
      cout << "x=";
      for(auto const v : vals) cout << v;
      cout << " z=" << zv << " ";

  	  uint res = 0;
  	  for(auto const v : vals)
  	  	res ^= v;
      if((res^c) == zv)
      	cout << "Correct value";
      else{
      	cout << "Incorrect value";
      	return false;
      }
      cout << endl;
    }
    cout << endl;


  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return true;
}

bool testRXCstXOR(uint const nvar){
	cout << "Test RX-diff propagation through XOR with constant for " << nvar << " variables" << endl;

	try {

	for(uint s = 0; s < nvar; s++){

		//Precompute the full table to check later
		uint64_t bound = (1ULL << nvar);
		vector<vector<vector<uint8_t>>> tablecheck(bound, vector<vector<uint8_t>>(bound, vector<uint8_t>(bound,0)));
		for(uint64_t dx = 0; dx < bound; dx++){
			for(uint64_t k = 0; k < bound; k++){
				for(uint64_t x = 0; x < bound; x++){
					uint64_t dy = CSHL(x^k,s,nvar) ^ CSHL(x,s,nvar)^dx^k;
					tablecheck[dx][k][dy] = 1;
				}
			}
		}
		uint64_t ctrCheck = 0;
		for(uint64_t dx = 0; dx < bound; dx++){
			for(uint64_t k = 0; k < bound; k++){
				for(uint64_t dy = 0; dy < bound; dy++){
					if(tablecheck[dx][k][dy] == 1)
						ctrCheck++;
				}
			}
		}

	    // Create an environment
	    GRBEnv env = GRBEnv();
	    env.set(GRB_IntParam_OutputFlag, 0);

	    // Create an empty model
	    GRBModel model = GRBModel(env);

	    // Create variables
	    vector<GRBVar> var(nvar);
	    vector<GRBVar> cst(nvar);
	    vector<GRBVar> res(nvar);
	    for(uint i = 0; i < nvar; i++){
	    	var[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i));
	    	cst[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "k"+to_string(i));
	    	res[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i));

	    }
	   	
	   	addRXCstXORConstr(model,var,cst,res,s);

	    model.set(GRB_IntParam_PoolSolutions, (1ULL << (3*nvar)));
	    model.set(GRB_IntParam_PoolSearchMode, 2);
	    model.optimize();

	    // Print number of solutions stored
	    uint nSolutions = model.get(GRB_IntAttr_SolCount);
	    cout << "Number of solutions found: " << nSolutions;
	    if(nSolutions == ctrCheck)
	    	cout << " : correct number of solutions" << endl;
	    else{
	    	cout << " : incorrect number of solutions, expected " << ctrCheck << endl;
	    	return false;
	    }

	    for (uint e = 0; e < nSolutions; e++) {
	      model.set(GRB_IntParam_SolutionNumber, e);

	      vector<uint> valx(nvar);
	      vector<uint> valk(nvar);
	      vector<uint> valy(nvar);
	      for(uint i = 0; i < nvar; i++){
	      	valx[i] = uint(round(var[i].get(GRB_DoubleAttr_Xn)));
	      	valk[i] = uint(round(cst[i].get(GRB_DoubleAttr_Xn)));
	      	valy[i] = uint(round(res[i].get(GRB_DoubleAttr_Xn)));
	      }
	      uint64_t xint = 0;
	      uint64_t kint = 0;
	      uint64_t yint = 0;
	      for(uint i = 0; i < nvar; i++){
	      	if(valx[i] != 0) xint |= (1ULL << i);
	      	if(valk[i] != 0) kint |= (1ULL << i);
	      	if(valy[i] != 0) yint |= (1ULL << i);
	      }
	      if(tablecheck[xint][kint][yint] == 0){
	      	cout << "Check failed, found solution x = " << xint << " k = " << kint << " y = " << yint << " with shift " << s << " but invalid" << endl;
	      	return false;
	      }

	    }
	    cout << endl;
	}
	return true;


  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return false;
}

bool testRXCstXOR_nonCstVar(uint const nvar){

	try {
	uint64_t bound = (1ULL << nvar);
	for(uint s = 0; s < nvar; s++){
		for(uint64_t cst = 0; cst < bound; cst++){

			//Write the constant in binary
			vector<uint> bincst(nvar);
			for(uint i = 0; i < nvar; i++)
				bincst[i] = (cst >> i)&1;

			//Precompute the full table to check later
			vector<vector<uint8_t>> tablecheck(bound, vector<uint8_t>(bound,0));
			for(uint64_t dx = 0; dx < bound; dx++){
				for(uint64_t x = 0; x < bound; x++){
					uint64_t dy = CSHL(x^cst,s,nvar) ^ CSHL(x,s,nvar)^dx^cst;
					tablecheck[dx][dy] = 1;
				}
			}
			uint64_t ctrCheck = 0;
			for(uint64_t dx = 0; dx < bound; dx++){
				for(uint64_t dy = 0; dy < bound; dy++){
					if(tablecheck[dx][dy] == 1)
						ctrCheck++;
				}
			}

		    // Create an environment
		    GRBEnv env = GRBEnv();
		    env.set(GRB_IntParam_OutputFlag, 0);

		    // Create an empty model
		    GRBModel model = GRBModel(env);

		    // Create variables
		    vector<GRBVar> var(nvar);
		    vector<GRBVar> res(nvar);
		    for(uint i = 0; i < nvar; i++){
		    	var[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i));
		    	res[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i));
		    }
		   	
		   	addRXCstXORConstr(model,var,bincst,res,s);

		    model.set(GRB_IntParam_PoolSolutions, (1ULL << (2*nvar)));
		    model.set(GRB_IntParam_PoolSearchMode, 2);
		    model.optimize();

		    // Print number of solutions stored
		    uint nSolutions = model.get(GRB_IntAttr_SolCount);
		    cout << "Number of solutions found: " << nSolutions;
		    if(nSolutions == ctrCheck)
		    	cout << " : correct number of solutions" << endl;
		    else{
		    	cout << " : incorrect number of solutions, expected " << ctrCheck << endl;
		    	return false;
		    }

		    for (uint e = 0; e < nSolutions; e++) {
				model.set(GRB_IntParam_SolutionNumber, e);

				vector<uint> valx(nvar);
				vector<uint> valy(nvar);
				for(uint i = 0; i < nvar; i++){
					valx[i] = uint(round(var[i].get(GRB_DoubleAttr_Xn)));
					valy[i] = uint(round(res[i].get(GRB_DoubleAttr_Xn)));
				}
				uint64_t xint = 0;
				uint64_t yint = 0;
				for(uint i = 0; i < nvar; i++){
					if(valx[i] != 0) xint |= (1ULL << i);
					if(valy[i] != 0) yint |= (1ULL << i);
				}
				if(tablecheck[xint][yint] == 0){
					cout << "Check failed, found solution x = " << xint << " k = " << cst << " y = " << yint << " with shift " << s << " but invalid" << endl;
					return false;
				}

		    }
		    cout << endl;
		}
	}

	return true;
  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return false;

}

bool testModAddValue(uint const n){

	uint mask = (1 << n) - 1;

	//Modelize and solve
	try {

    // Create an environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    vector<GRBVar> x(n);
    vector<GRBVar> y(n);
    vector<GRBVar> z(n);
    for(uint i = 0; i < n; i++){
    	x[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i));
    	y[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i));
    	z[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z"+to_string(i));
    }

    //Add mod add value constraints
    addModAddValueConstr(model,x,y,z,"c");

    //Optimize and find all solutions (should have 2**n total
    model.set(GRB_IntParam_PoolSolutions, (1ULL << (3*n)));
    model.set(GRB_IntParam_PoolSearchMode, 2);
    model.optimize();

    // Print number of solutions stored
    uint nSolutions = model.get(GRB_IntAttr_SolCount);
    cout << "Number of solutions found: " << nSolutions;
    if(nSolutions == (1ULL << (2*n)))
    	cout << " : correct number of solutions" << endl;
    else{
    	cout << " : incorrect number of solutions, expected " << (1ULL << (2*n)) << endl;
    	return false;
    }

    for (uint e = 0; e < nSolutions; e++) {
      model.set(GRB_IntParam_SolutionNumber, e);

      vector<uint> binx(n);
      vector<uint> biny(n);
      vector<uint> binz(n);
      for(uint i = 0; i < n; i++){
      	binx[i] = uint(round(x[i].get(GRB_DoubleAttr_Xn)));
      	biny[i] = uint(round(y[i].get(GRB_DoubleAttr_Xn)));
      	binz[i] = uint(round(z[i].get(GRB_DoubleAttr_Xn)));
      }

      uint valx = 0;
      uint valy = 0;
      uint valz = 0;
      for(uint i = 0; i < n; i++){
      	valx |= (binx[i] << i);
      	valy |= (biny[i] << i);
      	valz |= (binz[i] << i);
      }

      cout << "Sol ";
      if(e < 10) cout << " ";
      cout << e << " : ";
      cout << "x=" << valx << " y=" << valy << " z=" << valz;

      uint res = (valx + valy)&mask;
      if(res == valz) cout << " Correct value" << endl;
      else{
      	cout << " Incorrect value" << endl;
      	return false;
      }
    }
    cout << endl;
    return true;


  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return false;
}

bool testModAddConstValue(uint const n){

	uint mask = (1 << n) - 1;

	for(uint64_t cst = 0; cst < (1ULL << n); cst++){
		//Modelize and solve
		try {

	    // Create an environment
	    GRBEnv env = GRBEnv();
	    env.set(GRB_IntParam_OutputFlag, 0);

	    // Create an empty model
	    GRBModel model = GRBModel(env);

	    // Create variables
	    vector<GRBVar> x(n);
	    vector<GRBVar> z(n);
	    for(uint i = 0; i < n; i++){
	    	x[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i));
	    	z[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z"+to_string(i));
	    }

	    //Add mod add value constraints
	    addModAddConstValueConstr(model,x,cst,z,"c");

	    //Optimize and find all solutions 
	    model.set(GRB_IntParam_PoolSolutions, (1ULL << (2*n)));
	    model.set(GRB_IntParam_PoolSearchMode, 2);
	    model.optimize();

	    // Print number of solutions stored
	    uint nSolutions = model.get(GRB_IntAttr_SolCount);
	    cout << "Number of solutions found: " << nSolutions;
	    if(nSolutions == (1ULL << n))
	    	cout << " : correct number of solutions" << endl;
	    else{
	    	cout << " : incorrect number of solutions, expected " << (1ULL << n) << endl;
	    	return false;
	    }

	    for (uint e = 0; e < nSolutions; e++) {
	      model.set(GRB_IntParam_SolutionNumber, e);

	      vector<uint> binx(n);
	      vector<uint> binz(n);
	      for(uint i = 0; i < n; i++){
	      	binx[i] = uint(round(x[i].get(GRB_DoubleAttr_Xn)));
	      	binz[i] = uint(round(z[i].get(GRB_DoubleAttr_Xn)));
	      }

	      uint valx = 0;
	      uint valz = 0;
	      for(uint i = 0; i < n; i++){
	      	valx |= (binx[i] << i);
	      	valz |= (binz[i] << i);
	      }

	      cout << "Sol ";
	      if(e < 10) cout << " ";
	      cout << e << " : ";
	      cout << "x=" << valx << " y=" << cst << " z=" << valz;

	      uint res = (valx + cst)&mask;
	      if(res == valz) cout << " Correct value" << endl;
	      else{
	      	cout << " Incorrect value" << endl;
	      	return false;
	      }
	    }
	    cout << endl;
	   
	  } catch(GRBException e) {
	    cout << "Error code = " << e.getErrorCode() << endl;
	    cout << e.getMessage() << endl;
	    return false;
	  } catch(...) {
	    cout << "Exception during optimization" << endl;
	    return false;
	  }
	}
	return true;

}

bool testSpeckKS(uint const n){
	//Test vectors with full unrolling of the KS to verify each step
	//Tests done with Speck32/64, Speck64/128 and Speck128/256 only for now

	vector<uint64_t> key;
	vector<uint64_t> fullvalues;
	uint nbRound = 0;
	uint alpha = 8;
	uint beta = 3;
	uint m = 4;

	if(n == 16){
		key = vector<uint64_t>({0x0100,0x0908,0x1110,0x1918});
		fullvalues = vector<uint64_t>({0x100,0x1512,0x617d,0x1458,0x6919,0x77e2,0xc89,0xccdb,0xefea,0x4e33,0x76f4,0x5976,0xee8b,0xdb04,0x4617,0xf37e,0x87b4,0x8eca,0xed9b,0x3a52,0x8229,0xed64});
		alpha = 7;
		beta = 2;
		nbRound = fullvalues.size();
		m = key.size();
	}
	if(n == 24){
		key = vector<uint64_t>({0x020100,0x0a0908,0x121110});
		fullvalues = vector<uint64_t>({0x20100,0x1a0309,0xfa0d53,0xd37dc3,0x7549c5,0x7b02f3,0x8ee604,0x108772,0x2be5f4,0xdd6202,0xa901ff,0x48bbb,0xc9901f,0x435c51,0x11a899,0x4dfcb3,0xee268,0xfdb070,0x1d92e1,0xa28a9e,0x1d2e49,0x5a0e3a});
		nbRound = fullvalues.size();
		m = key.size();
	}
	else if(n == 32){
		key = vector<uint64_t>({0x03020100,0x0b0a0908,0x13121110,0x1b1a1918});
		fullvalues = vector<uint64_t>({0x3020100,0x131d0309,0xbbd80d53,0xd334df3,0x7fa43565,0x67e6ce55,0xe98cb3d2,0xaac76cbd,0x7f5951c8,0x3fa82c2,0x313533ad,0xdff70882,0x9e487c93,0xa934b928,0xdd2edef5,0x8be6388d,0x1f706b89,0x2b87aaf8,0x12d76c17,0x6eaccd6c,0x6a1ab912,0x10bc6bca,0x6057dd32,0xd3c9b381,0xb347813d,0x8c113c35,0xfe6b523a});
		nbRound = fullvalues.size();
		m = key.size();
	}
	else if(n == 64){
		// key = vector<uint64_t>({0x0706050403020100ULL,0x0f0e0d0c0b0a0908ULL,0x1716151413121110ULL,0x1f1e1d1c1b1a1918ULL});
		// fullvalues = vector<uint64_t>({0x706050403020100ULL,0x37253b31171d0309ULL,0xfe1588ce93d80d52ULL,0xe698e09f31334dfeULL,0xdb60f14bcbd834fdULL,0x2dafa7c34cc2c2f8ULL,0xfbb8e2705e64a1dbULL,0xdb6f99e4e383eaefULL,0x291a8d359c8ab92dULL,0xb653abee296e282ULL,0x604236be5c109d7fULL,0xb62528f28e15d89cULL,0x10419dd1d0b25f29ULL,0xfd71e73b9c69fff6ULL,0x8ea922047f976e93ULL,0x2e039afd398cffbcULL,0x9c9fcfef22c1072cULL,0x25fa8973ed55e6c9ULL,0x69819861a6b4280cULL,0x7b62d87498038f77ULL,0xf2351ece62e296feULL,0xa6d382d176ba05ffULL,0x8d96e66745b78726ULL,0xbe77397e9de6bf31ULL,0x35177f07af7d9479ULL,0xb86971c5e7815ff0ULL,0x7d77bfff103b45eaULL,0x9983914c82a1a11eULL,0x1e88e9b26e3307f5ULL,0x7a0068774fc7061bULL,0x1771e55c7df2b16fULL,0xa2cb5323bbf86418ULL,0x400303547ff5e38bULL,0xf4d26f589a56b276ULL});
		// nbRound = fullvalues.size();
		// m = key.size();

		key = vector<uint64_t>({0x0706050403020100ULL,0x0f0e0d0c0b0a0908ULL});
		fullvalues = vector<uint64_t>({0x706050403020100ULL,0x37253b31171d0309ULL,0xf91d89cc90c4085cULL,0xc6b1f07852cc7689ULL,0x14fcdf4f9c2d6f0ULL,0xb5fae1e4fe24cfd6ULL,0xa36d6954b0737cfeULL,0xf511691ea02f35f3ULL,0x5374abb75a2b455dULL,0x8dd5f6204ddcb2a5ULL,0xb243d7c9869cac18ULL,0x753e7a7c6660459eULL,0x78d648a3a5b0e63bULL,0x87152b23cbc0a8d2ULL,0xa8ff8b8c54a3b6f2ULL,0x4873be3c43b3ea79ULL,0x771ebffcbf05cb13ULL,0xe8a6bcaf25863d20ULL,0xe6c2ea8b5c520c93ULL,0x4d71b5c1ac5214f5ULL,0xdc60b2ae253070dcULL,0xb01d0abbe1fb9741ULL,0xd7987684a318b54aULL,0xa22c5282e600d319ULL,0xe029d67ebdf90048ULL,0x67559234c84efdbfULL,0x65173cf0cb01695cULL,0x24cf1f1879819519ULL,0x38a36ed2dbafb72aULL,0xded93cfe31bae304ULL,0xc53d18b91770b265ULL,0x2199c870db8ec93fULL});
		nbRound = fullvalues.size();
		m = key.size();
	}

	//Modelize and solve
	try {

    // Create an environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

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

	model.update();

	//Key schedule constraints
	addSpeckKSConstr(model,mk,k,alpha,beta);

	//Fix the master key to the test vector key value
	for(uint i = 0; i < m; i++){
		for(uint j = 0; j < n; j++){
			if((key[i] >> j)&1)
				model.addConstr(mk[i][j] == 1);
			else
				model.addConstr(mk[i][j] == 0);
		}
	}

	model.optimize();

    uint nSolutions = model.get(GRB_IntAttr_SolCount);
    if(nSolutions == 0){
    	cout << "Found no solution...." << endl;
    	return false;
    }
    else{
    	vector<uint64_t> keyval(nbRound);
    	for(uint i = 0; i < nbRound; i++){
    		if(i < 10) cout << " ";
    		cout << "k" << i << " : ";
    		for(uint j = 0; j < n; j++){
				uint v = uint(round(k[i][j].get(GRB_DoubleAttr_X)));
				cout << v;
				if(v == 1)
					keyval[i] |= (1ULL << j);
    		}
    		cout << " -> " << hex << setw(n/4) << setfill('0') << keyval[i] << dec;
    		if(keyval[i] == fullvalues[i])
    			cout << " correct value" << endl;
    		else{
    			cout << " incorrect value" << endl;
    			return false;
    		}
    	}
    }
    return true;

	} catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return false;
}

bool testModAddRXDiff(uint const n,
					  uint const gamma){

	uint64_t bound = (1ULL << n);
	uint64_t mask = bound-1;
	//Precompute all possibilities
	vector<vector<vector<uint64_t>>> RXDDT(bound, vector<vector<uint64_t>>(bound, vector<uint64_t>(bound,0)));

	for(uint64_t x = 0; x < bound; x++){
	 for(uint64_t y = 0; y < bound; y++){
	  uint64_t Rxpy = CSHL((x+y)&mask, gamma, n);
	  for(uint64_t a = 0; a < bound; a++){
	   uint64_t xRXa = CSHL(x,gamma,n) ^ a;
	   for(uint64_t b = 0; b < bound; b++){
	   	uint64_t d = ((xRXa + (CSHL(y,gamma,n) ^ b))&mask) ^ Rxpy;
	   	RXDDT[a][b][d]++;
	   }
	  }
	 }
	}

	//Count the number of valid transitions
	uint64_t validTransitions = 0;
	for(uint64_t a = 0; a < bound; a++){
	 for(uint64_t b = 0; b < bound; b++){
	  for(uint64_t d = 0; d < bound; d++){
	  	if(RXDDT[a][b][d] != 0)
			validTransitions++;
	  }
	 }
	}

	//Modelize and solve
	try {

    // Create an environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    vector<GRBVar> a(n);
    vector<GRBVar> b(n);
    vector<GRBVar> d(n);
    for(uint i = 0; i < n; i++){
    	a[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "a"+to_string(i));
    	b[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "b"+to_string(i));
    	d[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "d"+to_string(i));
    }

    //Add mod add value constraints
    addModAddRXDiffConstr(model,a,b,d,gamma);

    //Optimize and find all solutions (should have 2**n total
    model.set(GRB_IntParam_PoolSolutions, (1ULL << (3*n)));
    model.set(GRB_IntParam_PoolSearchMode, 2);
    model.optimize();

    // Print number of solutions stored
    uint64_t nSolutions = model.get(GRB_IntAttr_SolCount);
    cout << "Number of solutions found: " << nSolutions;
    if(nSolutions == validTransitions)
    	cout << " : correct number of solutions" << endl;
    else{
    	cout << " : incorrect number of solutions, expected " << validTransitions << endl;
    	return false;
    }

    for (uint e = 0; e < nSolutions; e++) {
      model.set(GRB_IntParam_SolutionNumber, e);

      vector<uint> bina(n);
      vector<uint> binb(n);
      vector<uint> bind(n);
      for(uint i = 0; i < n; i++){
      	bina[i] = uint(round(a[i].get(GRB_DoubleAttr_Xn)));
      	binb[i] = uint(round(b[i].get(GRB_DoubleAttr_Xn)));
      	bind[i] = uint(round(d[i].get(GRB_DoubleAttr_Xn)));
      }

      uint vala = 0;
      uint valb = 0;
      uint vald = 0;
      for(uint i = 0; i < n; i++){
      	vala |= (bina[i] << i);
      	valb |= (binb[i] << i);
      	vald |= (bind[i] << i);
      }

      cout << "Sol ";
      if(e < 10) cout << " ";
      cout << e << " : ";
      cout << "a=" << vala << " b=" << valb << " d=" << vald;

      if(RXDDT[vala][valb][vald] != 0) cout << " Correct value" << endl;
      else{
      	cout << " Incorrect value" << endl;
      	return false;
      }
    }
    cout << endl;
    return true;


  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return false;
}

bool testSpeckValue(uint const n,
					uint const m){
//Test the evaluation of Speck in value with MILP modeling
//Note: values for plaintext/ciphertext might have halves swapped compared to test vectors

	map<uint,map<uint,uint>> roundMap = {
		{16, {{4,22}}},
		{24, {{3,22}, {4,23}}},
		{32, {{3,26}, {4,27}}},
		{48, {{2,28}, {3,29}}},
		{64, {{2,32}, {3,33}, {4,34}}}
	};
	uint nbRound = roundMap[n][m];

	uint alpha = 8;
	uint beta = 3;
	if(n == 16){
		alpha = 7;
		beta = 2;
	}

	uint64_t px = 0;
	uint64_t py = 0;
	uint64_t cx = 0;
	uint64_t cy = 0;
	vector<uint64_t> pk;

	if(n == 16 && m == 4){
		px = 0x6574;
		py = 0x694c;
		cx = 0xa868;
		cy = 0x42f2;
		pk = vector<uint64_t>({0x0100,0x0908,0x1110,0x1918});
	}
	else if(n == 24 && m == 3){
		px = 0x20796c;
		py = 0x6c6172;
		cx = 0xc049a5;
		cy = 0x385adc;
		pk = vector<uint64_t>({0x020100,0x0a0908,0x121110});
	}
	else if(n == 64 && m == 2){
		px = 0x6c61766975716520;
		py = 0x7469206564616d20;
		cx = 0xa65d985179783265;
		cy = 0x7860fedf5c570d18;
		pk = vector<uint64_t>({0x0706050403020100ULL,0x0f0e0d0c0b0a0908ULL});
	}
	else{
		cout << "testSpeckValue not implemented with n = " << n << " m = " << m << endl;
		return false;
	}

	// Create an empty model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
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
	vector<vector<GRBVar>> z(nbRound, vector<GRBVar>(n));
	vector<vector<GRBVar>> y(nbRound+1, vector<GRBVar>(n));
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

	for(uint i = 0; i < nbRound; i++){
		vector<GRBVar> Saxi(n); //x[i] rotated to the right by alpha
		for(uint j = 0; j < n; j++)
			Saxi[j] = x[i][mod(j+alpha,n)];

		vector<GRBVar> Sbyi(n); //k[i] rotated to the left by beta
		for(uint j = 0; j < n; j++)
			Sbyi[j] = y[i][mod(j-beta,n)];

		//zi = Saxi + yi
		addModAddValueConstr(model,Saxi,y[i],z[i],"carry"+to_string(i));

		//x[i+1] = z[i] ^ k[i]
		//y[i+1] = Sbyi ^ x[i+1]
		for(uint j = 0; j < n; j++){
			addXORConstr(model, z[i][j], k[i][j], x[i+1][j]);
			addXORConstr(model, Sbyi[j], x[i+1][j], y[i+1][j]);
		}
	}

	//Add constraints for plaintext/key values
	for(uint i = 0; i < n; i++){
		if(((px >> i)&1) == 0)
			model.addConstr(x[0][i] == 0);
		else
			model.addConstr(x[0][i] == 1);

		if(((py >> i)&1) == 0)
			model.addConstr(y[0][i] == 0);
		else
			model.addConstr(y[0][i] == 1);
	}
	for(uint i = 0; i < m; i++){
		for(uint j = 0; j < n; j++){
			if(((pk[i] >> j)&1) == 0)
				model.addConstr(mk[i][j] == 0);
			else
				model.addConstr(mk[i][j] == 1);
		}
	}

	model.optimize();

	//Extract solution and check the ciphertext
	uint64_t nSolutions = model.get(GRB_IntAttr_SolCount);
	if(nSolutions == 0){
		cout << "No solution found" << endl;
		return false;
	}

	vector<uint64_t> bincx(n);
	vector<uint64_t> bincy(n);

	for(uint i = 0; i < n; i++){
		bincx[i] = uint(round(x[nbRound][i].get(GRB_DoubleAttr_X)));
		bincy[i] = uint(round(y[nbRound][i].get(GRB_DoubleAttr_X)));
	}

	uint64_t valcx = 0;
	uint64_t valcy = 0;
	for(uint i = 0; i < n; i++){
		valcx |= (bincx[i] << i);
		valcy |= (bincy[i] << i);
	}

	cout << " cx = " << hex << setw(n/4) << setfill('0') << valcx << dec;
	cout << " cy = " << hex << setw(n/4) << setfill('0') << valcy << dec;
	cout << endl;

	if(valcx == cx && valcy == cy)
		return true;
	else
		return false;
}

bool testLeaKS(uint const keySize){
	//Using test vectors from https://en.wikipedia.org/wiki/LEA_(cipher)

	vector<uint64_t> key;
	vector<vector<uint64_t>> fullvalues;

	if(keySize == 128){
		key = vector<uint64_t>({0x3c2d1e0f,0x78695a4b,0xb4a59687,0xf0e1d2c3});
		fullvalues = vector<vector<uint64_t>>({
			{0x3a0fd4,0x2497010,0x194f7db1,0x2497010,0x90d0883,0x2497010},
			{0x11fdcbb1,0x9e98e0c8,0x18b570cf,0x9e98e0c8,0x9dc53a79,0x9e98e0c8},
			{0xf30f7bb5,0x6d6628db,0xb74e5dad,0x6d6628db,0xa65e46d0,0x6d6628db},
			{0x74120631,0xdac9bd17,0xcd1ecf34,0xdac9bd17,0x540f76f1,0xdac9bd17},
			{0x662147db,0xc637c47a,0x46518932,0xc637c47a,0x23269260,0xc637c47a},
			{0xe4dd5047,0xf694285e,0xe1c2951d,0xf694285e,0x8ca5242c,0xf694285e},
			{0xbaf8e5ca,0x3e936cd7,0xfc7e5b1,0x3e936cd7,0xf1c8fa8c,0x3e936cd7},
			{0x5522b80c,0xee22ca78,0x8a6fa8b3,0xee22ca78,0x65637b74,0xee22ca78},
			{0x8a19279e,0x6fb40ffe,0x85c5f092,0x6fb40ffe,0x92cc9f25,0x6fb40ffe},
			{0x9dde584c,0xcb00c87f,0x4780ad66,0xcb00c87f,0xe61b5dcb,0xcb00c87f},
			{0x4fa10466,0xf728e276,0xd255411b,0xf728e276,0x656839ad,0xf728e276},
			{0x9250d058,0x51bd501f,0x1cb40dae,0x51bd501f,0x1abf218d,0x51bd501f},
			{0x21dd192d,0x77c644e2,0xcabfaa45,0x77c644e2,0x681c207d,0x77c644e2},
			{0xde7ac372,0x9436afd0,0x10331d80,0x9436afd0,0xf326fe98,0x9436afd0},
			{0xfb3ac3d4,0x93df660e,0x2f65d8a3,0x93df660e,0xdf92e761,0x93df660e},
			{0x27620087,0x265ef76e,0x4fb29864,0x265ef76e,0x2656ed1a,0x265ef76e},
			{0x227b88ec,0xd0b3fa6f,0xc86a08fd,0xd0b3fa6f,0xa864cba9,0xd0b3fa6f},
			{0xf1002361,0xe5e85fc3,0x1f0b0408,0xe5e85fc3,0x488e7ac4,0xe5e85fc3},
			{0xc65415d5,0x51e176b6,0xeca88bf9,0x51e176b6,0xedb89ece,0x51e176b6},
			{0x9b6fb99c,0x548254b,0x8de9f7c2,0x548254b,0xb6b4d146,0x548254b},
			{0x7257f134,0x6051a42,0x36bcef01,0x6051a42,0xb649d524,0x6051a42},
			{0xa540fb03,0x34b196e6,0xf7c80dad,0x34b196e6,0x71bc7dc4,0x34b196e6},
			{0x8fbee745,0xcf744123,0x907c0a60,0xcf744123,0x8215ec35,0xcf744123},
			{0xbf6adba,0xdf69029d,0x5b72305a,0xdf69029d,0xcb47c19f,0xdf69029d}});
	}
	else if(keySize == 192){
		key = vector<uint64_t>({0x3c2d1e0f,0x78695a4b,0xb4a59687,0xf0e1d2c3,0xc3d2e1f0,0x8796a5b4});
		fullvalues = vector<vector<uint64_t>>({
			{0x3a0fd4,0x2497010,0x194f7db1,0x90d0883,0x2ff5805a,0xc2580b27},
			{0x11fdcbb1,0x9e98e0c8,0x18b570cf,0x9dc53a79,0x5c145788,0x9771b5e5},
			{0xf30f7bb5,0x6d6628db,0xb74e5dad,0xa65e46d0,0x6f44da96,0xf643115f},
			{0x74120631,0xdac9bd17,0xcd1ecf34,0x540f76f1,0xaa1a5bdb,0xfbafaae7},
			{0x13f8a031,0x34f28728,0x31fdb409,0xe31481b,0xdf498117,0xcf9371f1},
			{0x967c312,0xb3484ec8,0x3aae5b3d,0x5a9714a0,0xb2d4dd5f,0x3a1fcdf7},
			{0xac47404,0x59e9e54d,0xa60dc00a,0x566139d3,0x898dce4f,0x582d72dd},
			{0x77f3ea4c,0xe2a73c8d,0xb8f1249a,0x6a172700,0xbc0e539c,0x2e46fdbb},
			{0xb4e0e98a,0x3d028c05,0xb8d3a050,0xdbd67bef,0xdf675c7a,0x99eefbb0},
			{0xe68584f6,0xce31ef45,0x96c105ac,0x2a1be677,0x9d72b8b0,0x33cecc54},
			{0xc22ffd76,0x1ab7167e,0x42bb3060,0x7da517f5,0x4aa0e8d3,0xa070c3c},
			{0xe200a765,0xc2be17b3,0x7f22543f,0x3e4eb7a1,0xc992a6f4,0xa783c823},
			{0xc13cc747,0xffcc8185,0x66514e9e,0xe4ccc199,0xcd5c766d,0xa004f676},
			{0x1d3a1fa6,0xd46894ec,0xf49c33e6,0x782fda7e,0x1fe6346c,0xffe981c},
			{0x78b97c3d,0x956e8ee8,0x49ab721c,0x2672138a,0x37ea242,0xce5fe8a4},
			{0x225f7158,0x32d83e3e,0xe118f6aa,0x1fb83751,0x4d27715c,0xed2fba4e},
			{0x8dfbc56d,0xe0a907db,0xe4af091c,0x5e123225,0xd0e8d2e1,0xcc4501fb},
			{0x8422a8f0,0x46a12f92,0x415152ad,0xf55417f5,0x38738248,0xc6e29ded},
			{0x5723715e,0xabfa788c,0xc3646af7,0x64af9186,0x8fc855ec,0x2bc36989},
			{0x5e6b28e3,0xe0f5f592,0xeb3dd108,0x551012a,0x50e4221d,0x97e85c0f},
			{0x4e258e14,0x92298f0b,0x771269c3,0x6f934254,0xc0933b6b,0x421159b8},
			{0xd76953f4,0x6a3e36be,0x53b656fb,0x610c22e0,0x9f399330,0xacf7e7e9},
			{0xfe0b573b,0xcbb73085,0x89ed67fc,0x77014cef,0xe1b8431f,0xba1b4105},
			{0x6de3450,0xb3f5b2fe,0xdf1cec27,0xfb22bd10,0x8e3de6fe,0x3d4acd27},
			{0xc5444873,0x5bec968b,0x8b2af393,0x11e2f6ca,0x9cb3694f,0x94c56b91},
			{0x939a1a93,0x27f101bb,0x5381bae7,0x48ebd1b1,0xf6d5fca7,0xca24bbc},
			{0x7b03490b,0xde00acfb,0xc7f8abfe,0x410a14c1,0xd37932a9,0x14029327},
			{0xbd948525,0x2c75004d,0xc52486d5,0xf07e2fa,0x1963e1fd,0x882719c3}});
	}
	else if(keySize == 256){
		key = vector<uint64_t>({0x3c2d1e0f,0x78695a4b,0xb4a59687,0xf0e1d2c3,0xc3d2e1f0,0x8796a5b4,0x4b5a6978,0xf1e2d3c});
		fullvalues = vector<vector<uint64_t>>({
			{0x3a0fd4,0x2497010,0x194f7db1,0x90d0883,0x2ff5805a,0xc2580b27},
			{0xa83e7ef9,0x53eca29,0xd359f988,0x8101a243,0x9bbf34b3,0x9228434f},
			{0x2efee506,0x8b5f7bd4,0x9991e811,0x72dbc20c,0x2384c97f,0xcefee47f},
			{0xc571782c,0xda90b1,0xb940a552,0x5db79619,0x4bc9a125,0x5d08a419},
			{0x72de26cc,0xd69bc26f,0x46a7f207,0x66ff4d81,0xa87862fc,0xa5f63601},
			{0x7909c4fa,0xf3f93651,0x72cb0bcd,0xae69b2e3,0x80f2ca4b,0xf13efcce},
			{0x7869db69,0x6b7a5b8e,0xfefbf6b1,0xec608c8e,0x76e9d5d2,0x13ca4bf6},
			{0xc5eeec7a,0xaa42a59d,0x1f22cd00,0xfdd92bdc,0xd6bbe3e8,0x15d459ec},
			{0xcda7632a,0x9cf01bef,0x6596e261,0x8c1de14c,0x1127c3b8,0x48b3f629},
			{0x3723d0e1,0xfc0317ec,0x3fdd5378,0x201ae1d,0xe55db65e,0xe4c84dbc},
			{0x3633db3f,0xe4c24fc2,0xbb1e1fd7,0xa339425c,0xfe3e1bdf,0xd61c808d},
			{0xbdca3449,0xbeb8aa4e,0x145a9687,0xeb6fcd87,0x8b88ca72,0x7677a84b},
			{0xd11005e9,0x558275c5,0xbc742819,0x3f17e888,0x20fcb71f,0x60886959},
			{0x8d9446c4,0x67d2d167,0x855a6aef,0x69ea517c,0x36e48e11,0xd3f4e86},
			{0xbb0ede65,0xcceecc06,0xefc9c49f,0x44902261,0xbd8549c0,0xa7e7f682},
			{0x772101e6,0xb4b9a250,0x6faa7b73,0x7318b792,0x1e57e751,0xfd43b41c},
			{0x4ec21b5f,0xdcfbf30b,0xa4046947,0xbe0e781c,0xd74e21ac,0x6b1f5d22},
			{0xe8b8e02b,0x4a662d2d,0xb50f9ca9,0x1c98c69,0x9eb28089,0x216cfd3f},
			{0x92f0126b,0x7b9961aa,0x581f94ac,0xab4be6dd,0xc2a91af5,0xfb4e8e0c},
			{0x4c2c8f04,0x81a45991,0x1fcb946c,0xbccbb5b5,0x808899cb,0x8c1b2f89},
			{0x192061be,0x78e5cf04,0xf239ab5c,0xe8471e86,0x9e6217c7,0xe5fdf35c},
			{0x83c3150d,0x766887f8,0xa1092ac7,0x6aa6f41d,0x16e200f9,0x6bdc26ca},
			{0x52345706,0xdb70d6af,0xa8d8ffeb,0x492ee661,0x4cd1e991,0xd75d8352},
			{0x85a9c5fb,0x1e0f569e,0x7ff7c600,0x3f36a1d8,0xe406ad00,0x4ded8f16},
			{0x512bb2f4,0x772b192c,0x2e6168bd,0x76af67e1,0xd893a786,0x3e276f69},
			{0xd11ee3ad,0xb7f8c612,0xd3b19318,0x89fee4db,0xb6c3aedd,0x5420f90},
			{0x4f662f0,0x8fb41a6c,0x2f42dd5e,0xa8ad1839,0x46474e45,0x46418de0},
			{0x351550c8,0x668014f6,0x4924365,0x5f353d6f,0x4eba8d76,0x924a4318},
			{0x5aba711c,0xa36b1398,0x5b3e7bf4,0x7b3a2cf9,0x1d006ebe,0xd5683e5},
			{0x4f56916f,0x215dccd2,0x9f57886f,0x876d1357,0x46013d49,0x2a4932a3},
			{0xaa285691,0xebefe7d3,0xe960e64b,0xdd893f0f,0x6a234412,0x495d13c9},
			{0x71c683e8,0x8069dfd0,0x6c1a501d,0x699418,0x262142f0,0xa91a7393}});
	}

	uint nbRound = fullvalues.size();
	uint m = key.size();
	uint n = 32;

		//Modelize and solve
	try {

    // Create an environment
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

	// Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables
    //Key variables
    vector<vector<GRBVar>> mk(m, vector<GRBVar>(32)); //Special variables for the master key, for easier use
    vector<vector<vector<GRBVar>>> k(nbRound, vector<vector<GRBVar>>(6, vector<GRBVar>(32)));
	for(uint i = 0; i < nbRound; i++){
		for(uint l = 0; l < 6; l++){
			for(uint j = 0; j < n; j++)
				k[i][l][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "k"+to_string(i)+"_"+to_string(l)+"_"+to_string(j));
		}
	}
	for(uint i = 0; i < m; i++){
		for(uint j = 0; j < n; j++)
			mk[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "mk"+to_string(i)+"_"+to_string(j));
	}

	model.update();

	//Key schedule constraints
	addLeaKSConstr(model,mk,k);

	//Fix the master key to the test vector key value
	for(uint i = 0; i < m; i++){
		for(uint j = 0; j < 32; j++){
			if((key[i] >> j)&1)
				model.addConstr(mk[i][j] == 1);
			else
				model.addConstr(mk[i][j] == 0);
		}
	}

	model.optimize();

	uint nSolutions = model.get(GRB_IntAttr_SolCount);
    if(nSolutions == 0){
    	cout << "Found no solution...." << endl;
    	return false;
    }
    else{
    	vector<vector<uint64_t>> keyval(nbRound, vector<uint64_t>(6,0));
    	for(uint r = 0; r < nbRound; r++){
    		if(r < 10) cout << " ";
    		cout << "k" << r << " : " << endl;
    		for(uint i = 0; i < 6; i++){
    			for(uint j = 0; j < 32; j++){
    				uint v = uint(round(k[r][i][j].get(GRB_DoubleAttr_X)));
    				cout << v;
    				if(v == 1)
    					keyval[r][i] |= (1ULL << j);
    			}
    			cout << " -> " << hex << setw(n/4) << setfill('0') << keyval[r][i] << dec;
    			if(keyval[r][i] == fullvalues[r][i])
    				cout << " correct value" << endl;
    			else{
    				cout << " incorrect value" << endl;
    				return false;
    			}
    		}
    	}
    }
    return true;



	} catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return false;
}

bool testLeaValue(uint const keySize){

	uint nbRound = 0;
	vector<uint64_t> p;
	vector<uint64_t> c;
	vector<uint64_t> key;


	if(keySize == 128){
		nbRound = 24;
		p = vector<uint64_t>({0x13121110,0x17161514,0x1b1a1918,0x1f1e1d1c});
		c = vector<uint64_t>({0x354ec89f,0x18c6c628,0xa7c73255,0xfd8b6404});
		key = vector<uint64_t>({0x3c2d1e0f,0x78695a4b,0xb4a59687,0xf0e1d2c3});
	}
	else if(keySize == 192){
		nbRound = 28;
		p = vector<uint64_t>({0x23222120,0x27262524,0x2b2a2928,0x2f2e2d2c});
		c = vector<uint64_t>({0x325eb96f,0x871bad5a,0x35f5dc8c,0xf2c67476});
		key = vector<uint64_t>({0x3c2d1e0f,0x78695a4b,0xb4a59687,0xf0e1d2c3,0xc3d2e1f0,0x8796a5b4});
	}
	else if(keySize == 256){
		nbRound = 32;
		p = vector<uint64_t>({0x33323130,0x37363534,0x3b3a3938,0x3f3e3d3c});
		c = vector<uint64_t>({0xf6af51d6,0xc189b147,0xca00893a,0x97e1f927});
		key = vector<uint64_t>({0x3c2d1e0f,0x78695a4b,0xb4a59687,0xf0e1d2c3,0xc3d2e1f0,0x8796a5b4,0x4b5a6978,0xf1e2d3c});
	}
	


	try{
	// Create an empty model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
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

    model.update();

    //Key schedule constraints
    addLeaKSConstr(model,mk,k);

    //Round functions
    for(uint r = 0; r < nbRound; r++){
    	//Step 1
    	for(uint i = 0; i < 32; i++){
	    	addXORConstr(model, x[r][0][i], k[r][0][i], y[r][0][i]);
	    	addXORConstr(model, x[r][1][i], k[r][1][i], y[r][1][i]);
	    	addXORConstr(model, x[r][1][i], k[r][2][i], y[r][2][i]);
	    	addXORConstr(model, x[r][2][i], k[r][3][i], y[r][3][i]);
	    	addXORConstr(model, x[r][2][i], k[r][4][i], y[r][4][i]);
	    	addXORConstr(model, x[r][3][i], k[r][5][i], y[r][5][i]);
	    }

    	//Step 2:
    	vector<GRBVar> out(32);
    	for(uint i = 0; i < 32; i++)
    		out[i] = x[r+1][0][(i+9)%32];
    	addModAddValueConstr(model,y[r][0],y[r][1],out,"carry0"+to_string(r));

    	for(uint i = 0; i < 32; i++)
    		out[i] = x[r+1][1][mod(i-5,32)];
    	addModAddValueConstr(model,y[r][2],y[r][3],out,"carry1"+to_string(r));

    	for(uint i = 0; i < 32; i++)
    		out[i] = x[r+1][2][mod(i-3,32)];
    	addModAddValueConstr(model,y[r][4],y[r][5],out,"carry2"+to_string(r));

    	for(uint i = 0; i < 32; i++)
    		model.addConstr(x[r+1][3][i] == x[r][0][i]);
    }

    //Add constraints for plaintext/key values
    for(uint i = 0; i < 4; i++){
    	for(uint j = 0; j < 32; j++){
    		if(((p[i] >> j)&1) == 0)
    			model.addConstr(x[0][i][j] == 0);
    		else
    			model.addConstr(x[0][i][j] == 1);
    	}
    }
    for(uint i = 0; i < m; i++){
    	for(uint j = 0; j < 32; j++){
    		if(((key[i] >> j)&1) == 0)
    			model.addConstr(mk[i][j] == 0);
    		else
    			model.addConstr(mk[i][j] == 1);
    	}
    }

    model.optimize();

    //Extract solution and check the ciphertext
	uint64_t nSolutions = model.get(GRB_IntAttr_SolCount);
	if(nSolutions == 0){
		cout << "No solution found" << endl;
		return false;
	}

	vector<vector<uint64_t>> binc(4,vector<uint64_t>(32));
	for(uint i = 0; i < 4; i++){
		for(uint j = 0; j < 32; j++)
			binc[i][j] = uint(round(x[nbRound][i][j].get(GRB_DoubleAttr_X)));
	}

	vector<uint64_t> valc(4,0);
	cout << "c = ";
	for(uint i = 0; i < 4; i++){
		for(uint j = 0; j < 32; j++)
			valc[i] |= (binc[i][j] << j);
		cout << hex << setw(32/4) << setfill('0') << valc[i] << dec << " ";	
		if(valc[i] != c[i]){
			cout << "Wrong ciphertext value" << endl;
			return false;
		}
	}
	return true;



	} catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }
  return false;

}

std::pair<uint32_t,uint32_t> alzetteCore(std::vector<uint> const & params,
										 uint32_t const cst,
										 uint32_t const inx,
										 uint32_t const iny){

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

bool testAlzetteValue(std::vector<uint> const & params,
					  uint const cst){

	//Pick a random input
	std::random_device rd;
	std::mt19937_64 prng(rd());
	std::uniform_int_distribution<uint32_t> randuint32;
	uint32_t inputx = randuint32(prng);
	uint32_t inputy = randuint32(prng);
	auto expectedOutput = alzetteCore(params,cst,inputx,inputy);

	uint nbRound = 4;
	//Write the constant in binary
	vector<uint> bincst(32);
	for(uint i = 0; i < 32; i++)
		bincst[i] = (cst >> i)&1;

	try{

	// Create an empty model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
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

    for(uint r = 0;  r < nbRound; r++){

    	vector<GRBVar> tmp(32);

    	//z[r] = x[r] + (y[r] >>> params[2r])
    	for(uint i = 0; i < 32; i++)
    		tmp[i] = y[r][(i+params[2*r])%32]; //tmp = y[r] >>> params[2r]
    	addModAddValueConstr(model, x[r], tmp, z[r], "carry"+to_string(r));

    	//x[r+1] = z[r] ^ cst
    	for(uint i = 0; i < 32; i++){
    		if(bincst[i] == 0)
    			model.addConstr(x[r+1][i] == z[r][i]);
    		else
    			model.addConstr(x[r+1][i] == (1 - z[r][i]));
    	}

    	//y[r+1] = y[r] ^ (z[r] >>> params[2r+1])
    	for(uint i = 0; i < 32; i++)
    		tmp[i] = z[r][(i+params[2*r+1])%32]; //tmp = z[r] >>> params[2r+1]
    	for(uint i = 0; i < 32; i++)
    		addXORConstr(model, y[r][i], tmp[i], y[r+1][i]);
    }

    //Add constraints on the input
    for(uint i = 0; i < 32; i++){
    	if(((inputx >> i)&1) == 0)
    		model.addConstr(x[0][i] == 0);
    	else
    		model.addConstr(x[0][i] == 1);

    	if(((inputy >> i)&1) == 0)
    		model.addConstr(y[0][i] == 0);
    	else
    		model.addConstr(y[0][i] == 1);
    }

    model.update();
    model.optimize();

    //Extract solution and check the ciphertext
	uint64_t nSolutions = model.get(GRB_IntAttr_SolCount);
	if(nSolutions == 0){
		cout << "No solution found" << endl;
		return false;
	}
	else{
		vector<uint> binx(32);
		vector<uint> biny(32);
		for(uint i = 0; i < 32; i++){
			binx[i] = uint(round(x[nbRound][i].get(GRB_DoubleAttr_X)));
			biny[i] = uint(round(y[nbRound][i].get(GRB_DoubleAttr_X)));
		}
		uint32_t resultx = 0;
		uint32_t resulty = 0;
		for(uint i = 0; i < 32; i++){
			resultx |= (binx[i] << i);
			resulty |= (biny[i] << i);
		}

		if(resultx != expectedOutput.first || resulty != expectedOutput.second){
			cout << "Wrong result, input x = ";
			cout << hex << setw(32/4) << setfill('0') << inputx << " y = ";	
			cout << hex << setw(32/4) << setfill('0') << inputy << endl;
			cout << "Expected output x = ";
			cout << hex << setw(32/4) << setfill('0') << expectedOutput.first << " y = ";
			cout << hex << setw(32/4) << setfill('0') << expectedOutput.second << endl;
			cout << "Obtained output x = ";
			cout << hex << setw(32/4) << setfill('0') << resultx << " y = ";
			cout << hex << setw(32/4) << setfill('0') << resulty << endl;
			cout << dec;
			return false;
		}
	}
	return true;


	} catch(GRBException e) {
cout << "Error code = " << e.getErrorCode() << endl;
cout << e.getMessage() << endl;
} catch(...) {
cout << "Exception during optimization" << endl;
}
return false;

}

bool testXteaKS(std::vector<uint32_t> const & key){

	uint nbRound = 64;
	uint keySize = 128;

	//Precompute the keyschedule with the given key
	vector<uint32_t> roundConstants({0x0,0x9e3779b9,0x9e3779b9,0x3c6ef372,0x3c6ef372,0xdaa66d2b,0xdaa66d2b,0x78dde6e4,0x78dde6e4,0x1715609d,0x1715609d,0xb54cda56,0xb54cda56,0x5384540f,0x5384540f,0xf1bbcdc8,0xf1bbcdc8,0x8ff34781,0x8ff34781,0x2e2ac13a,0x2e2ac13a,0xcc623af3,0xcc623af3,0x6a99b4ac,0x6a99b4ac,0x8d12e65,0x8d12e65,0xa708a81e,0xa708a81e,0x454021d7,0x454021d7,0xe3779b90,0xe3779b90,0x81af1549,0x81af1549,0x1fe68f02,0x1fe68f02,0xbe1e08bb,0xbe1e08bb,0x5c558274,0x5c558274,0xfa8cfc2d,0xfa8cfc2d,0x98c475e6,0x98c475e6,0x36fbef9f,0x36fbef9f,0xd5336958,0xd5336958,0x736ae311,0x736ae311,0x11a25cca,0x11a25cca,0xafd9d683,0xafd9d683,0x4e11503c,0x4e11503c,0xec48c9f5,0xec48c9f5,0x8a8043ae,0x8a8043ae,0x28b7bd67,0x28b7bd67,0xc6ef3720});
	vector<uint32_t> fullvalues(nbRound);
	for(uint i = 0; i < nbRound; i++){
		if(i%2 == 0)
			fullvalues[i] = roundConstants[i] + key[roundConstants[i]&3];
		else
			fullvalues[i] = roundConstants[i] + key[(roundConstants[i] >> 11)&3];
	}

	try{

	// Create an empty model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
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

    model.update();

    //Key schedule constraints
    addXteaKSConstr(model,mk,k);

    //Fix the key to the test vector key value
    for(uint i = 0; i < 4; i++){
    	for(uint j = 0; j < 32; j++){
    		if((key[i] >> j)&1)
				model.addConstr(mk[i][j] == 1);
			else
				model.addConstr(mk[i][j] == 0);
		}
    }

    model.optimize();

    uint nSolutions = model.get(GRB_IntAttr_SolCount);
    if(nSolutions == 0){
    	cout << "Found no solution...." << endl;
    	return false;
    }
    else{
    	vector<uint64_t> keyval(nbRound, 0);
    	for(uint r = 0; r < nbRound; r++){
    		if(r < 10) cout << " ";
    		cout << "k" << r << " : " << endl;
			for(uint j = 0; j < 32; j++){
				uint v = uint(round(k[r][j].get(GRB_DoubleAttr_X)));
				cout << v;
				if(v == 1)
					keyval[r] |= (1ULL << j);
			}
			cout << " -> " << hex << setw(8) << setfill('0') << keyval[r] << dec;
			if(keyval[r] == fullvalues[r])
				cout << " correct value" << endl;
			else{
				cout << " incorrect value" << endl;
				return false;
			}
    	}
    }
    return true;

	}  catch(GRBException e) {
cout << "Error code = " << e.getErrorCode() << endl;
cout << e.getMessage() << endl;
} catch(...) {
cout << "Exception during optimization" << endl;
}
return false;

}

/* take 64 bits of data in v[0] and v[1] and 128 bits of key[0] - key[3] 
Code taken from wikipedia, no official test vectors but it seems to match several other implementations found online*/
void xteaEncrypt(unsigned int num_rounds, uint32_t v[2], uint32_t const key[4]){
	unsigned int i;
    uint32_t v0=v[0], v1=v[1], sum=0, delta=0x9E3779B9;
    for (i=0; i < num_rounds; i++) {
        v0 += (((v1 << 4) ^ (v1 >> 5)) + v1) ^ (sum + key[sum & 3]);
        sum += delta;
        v1 += (((v0 << 4) ^ (v0 >> 5)) + v0) ^ (sum + key[(sum>>11) & 3]);
    }
    v[0]=v0; v[1]=v1;
}

bool testXteaValue(){

	uint nbRound = 64;
	uint keySize = 128;

	//Pick a random input and key
	std::random_device rd;
	std::mt19937_64 prng(rd());
	std::uniform_int_distribution<uint32_t> randuint32;
	uint32_t inputx = randuint32(prng);
	uint32_t inputy = randuint32(prng);
	uint32_t inputkey[4];
	for(uint i = 0; i < 4; i++)
		inputkey[i] = randuint32(prng);
	uint32_t expectedOutput[2] = {inputx,inputy};

	cout << "Expected input x = ";
	cout << hex << setw(32/4) << setfill('0') << expectedOutput[0] << " y = ";
	cout << hex << setw(32/4) << setfill('0') << expectedOutput[1] << endl;
	xteaEncrypt(32,expectedOutput,inputkey);
	cout << "Expected output x = ";
	cout << hex << setw(32/4) << setfill('0') << expectedOutput[0] << " y = ";
	cout << hex << setw(32/4) << setfill('0') << expectedOutput[1] << endl;
	cout << dec;

	try{

	// Create an empty model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
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
    vector<vector<GRBVar>> a(nbRound, vector<GRBVar>(32));
    vector<vector<GRBVar>> b(nbRound, vector<GRBVar>(32));
    vector<vector<GRBVar>> c(nbRound, vector<GRBVar>(32));
    for(uint i = 0; i < nbRound; i++){
    	for(uint j = 0; j < 32; j++){
    		a[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "a"+to_string(i)+"_"+to_string(j));
    		b[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "b"+to_string(i)+"_"+to_string(j));
    		c[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "c"+to_string(i)+"_"+to_string(j));
    	}
    }

    model.update();

    //Key schedule constraints
    addXteaKSConstr(model,mk,k);

    /*
	input x[i] y[i]
	a[i] = (y[i] << 4) ^ (y[i] >> 5)
	b[i] = a[i] + y[i]
	c[i] = b[i] ^ roundkey[i]
	y[i+1] = x[i] + c[i]
	x[i+1] = y[i]
	output x[i+1] y[i+1]
	*/

    //Round functions
    for(uint r = 0; r < nbRound; r++){

    	//Step1: a[i] = (y[i] << 4) ^ (y[i] >> 5)
    	//It's a shift, not a rotation, slightly annoying
    	//  y27 y26 y25 y24 y23 y22 ... y0  0  0  0  0
    	//^   0   0   0   0   0 y31 ... y9 y8 y7 y6 y5
    	//= a31 a30 a29 a28 a27 a26 ... a4 a3 a2 a1 a0
    	for(uint i = 0; i < 4; i++)
    		model.addConstr(a[r][i] == y[r][i+5]);
    	for(uint i = 4; i < 27; i++)
    		addXORConstr(model, y[r][i-4], y[r][i+5], a[r][i]);
    	for(uint i = 27; i < 32; i++)
    		model.addConstr(a[r][i] == y[r][i-4]);

    	//Step2: b[i] = a[i] + y[i]
    	addModAddValueConstr(model,a[r],y[r],b[r],"carryb"+to_string(r));

    	//Step3: c[i] = b[i] ^ roundkey[i]
    	for(uint i = 0; i < 32; i++)
    		addXORConstr(model, b[r][i], k[r][i], c[r][i]);

    	//Step4: y[i+1] = x[i] + c[i]
    	addModAddValueConstr(model,x[r],c[r],y[r+1],"carryc"+to_string(r));

    	//Step5: x[i+1] = y[i]
    	for(uint i = 0; i < 32; i++)
    		model.addConstr(x[r+1][i] == y[r][i]);
    }

    //Add constraints on the input
    for(uint i = 0; i < 32; i++){
    	if(((inputx >> i)&1) == 0)
    		model.addConstr(x[0][i] == 0);
    	else
    		model.addConstr(x[0][i] == 1);

    	if(((inputy >> i)&1) == 0)
    		model.addConstr(y[0][i] == 0);
    	else
    		model.addConstr(y[0][i] == 1);
    }

    //Add constraints on the key
    for(uint i = 0; i < 4; i++){
    	for(uint j = 0; j < 32; j++){
    		if(((inputkey[i] >> j)&1) == 0)
    			model.addConstr(mk[i][j] == 0);
    		else
    			model.addConstr(mk[i][j] == 1);
    	}
    }

    model.update();
    model.optimize();

    //Extract solution and check the ciphertext
	uint64_t nSolutions = model.get(GRB_IntAttr_SolCount);
	if(nSolutions == 0){
		cout << "No solution found" << endl;
		return false;
	}
	else{
		vector<uint> binx(32);
		vector<uint> biny(32);
		for(uint i = 0; i < 32; i++){
			binx[i] = uint(round(x[nbRound][i].get(GRB_DoubleAttr_X)));
			biny[i] = uint(round(y[nbRound][i].get(GRB_DoubleAttr_X)));
		}
		uint32_t resultx = 0;
		uint32_t resulty = 0;
		for(uint i = 0; i < 32; i++){
			resultx |= (binx[i] << i);
			resulty |= (biny[i] << i);
		}

		if(resultx != expectedOutput[0] || resulty != expectedOutput[1]){
			cout << "Wrong result, input x = ";
			cout << hex << setw(32/4) << setfill('0') << inputx << " y = ";	
			cout << hex << setw(32/4) << setfill('0') << inputy << endl;
			cout << "Expected output x = ";
			cout << hex << setw(32/4) << setfill('0') << expectedOutput[0] << " y = ";
			cout << hex << setw(32/4) << setfill('0') << expectedOutput[1] << endl;
			cout << "Obtained output x = ";
			cout << hex << setw(32/4) << setfill('0') << resultx << " y = ";
			cout << hex << setw(32/4) << setfill('0') << resulty << endl;
			cout << dec;
			return false;
		}
	}
	return true;

}  catch(GRBException e) {
cout << "Error code = " << e.getErrorCode() << endl;
cout << e.getMessage() << endl;
} catch(...) {
cout << "Exception during optimization" << endl;
}
return false;

}

bool testSHLRXDiff(int const n,
				   int const alpha,
				   int const gamma){

	// Create an empty model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);

	uint64_t bound = (1ULL << n);
	//precompute the full table
	vector<vector<bool>> table(bound, vector<bool>(bound,false));
	for(uint64_t x0 = 0; x0 < bound; x0++){
		uint64_t rx0 = CSHL(x0, gamma,n);
		uint64_t ry0 = CSHL(shl(x0,alpha,n), gamma,n);
		for(uint64_t x1 = 0; x1 < bound; x1++){
			uint64_t y1 = shl(x1,alpha,n);
			uint64_t deltaIn = rx0 ^ x1;
			uint64_t deltaOut = ry0 ^ y1;
			table[deltaIn][deltaOut] = true;
		}
	}
	uint64_t ctrCheck = 0;
	for(uint64_t deltaIn = 0; deltaIn < bound; deltaIn++){
		for(uint64_t deltaOut = 0; deltaOut < bound; deltaOut++){
			if(table[deltaIn][deltaOut])
				ctrCheck++;
		}
	}


	try{
    GRBModel model = GRBModel(env);

    //State variables
    vector<GRBVar> x(n);
    vector<GRBVar> y(n);
    for(int i = 0; i < n; i++){
		x[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i));
		y[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i));
    }
    model.update();

    //SHL constraints
	addSHLRXDiffConstr(model,x,y,alpha,gamma);

    model.update();
    model.optimize();

    model.set(GRB_IntParam_PoolSolutions, 1ULL << (2*n));
    model.set(GRB_IntParam_PoolSearchMode, 2);
    model.optimize();

    // Print number of solutions stored
    uint nSolutions = model.get(GRB_IntAttr_SolCount);
    cout << "Number of solutions found: " << nSolutions;
    if(nSolutions == ctrCheck)
    	cout << " : correct number of solutions" << endl;
    else{
    	cout << " : incorrect number of solutions, expected " << ctrCheck << endl;
    	return false;
    }

    for (uint e = 0; e < nSolutions; e++) {
      model.set(GRB_IntParam_SolutionNumber, e);

      vector<uint> valx(n);
      vector<uint> valy(n);
      for(int i = 0; i < n; i++){
      	valx[i] = uint(round(x[i].get(GRB_DoubleAttr_Xn)));
      	valy[i] = uint(round(y[i].get(GRB_DoubleAttr_Xn)));
      }
      uint64_t xint = 0;
      uint64_t yint = 0;
      for(int i = 0; i < n; i++){
      	if(valx[i] != 0) xint |= (1ULL << i);
      	if(valy[i] != 0) yint |= (1ULL << i);
      }
      if(!table[xint][yint]){
      	cout << "Check failed, found solution x = " << xint << " y = " << yint << " with alpha = " << alpha + " gamma = " << gamma << endl;
      	return false;
      }

    }
    cout << endl;
	return true;

}  catch(GRBException e) {
cout << "Error code = " << e.getErrorCode() << endl;
cout << e.getMessage() << endl;
} catch(...) {
cout << "Exception during optimization" << endl;
}
return false;
}

bool testSHRRXDiff(int const n,
				   int const beta,
				   int const gamma){

	// Create an empty model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);

	uint64_t bound = (1ULL << n);
	//precompute the full table
	vector<vector<bool>> table(bound, vector<bool>(bound,false));
	for(uint64_t x0 = 0; x0 < bound; x0++){
		uint64_t rx0 = CSHL(x0, gamma,n);
		uint64_t ry0 = CSHL(shr(x0,beta,n), gamma,n);
		for(uint64_t x1 = 0; x1 < bound; x1++){
			uint64_t y1 = shr(x1,beta,n);
			uint64_t deltaIn = rx0 ^ x1;
			uint64_t deltaOut = ry0 ^ y1;
			table[deltaIn][deltaOut] = true;
		}
	}
	uint64_t ctrCheck = 0;
	for(uint64_t deltaIn = 0; deltaIn < bound; deltaIn++){
		for(uint64_t deltaOut = 0; deltaOut < bound; deltaOut++){
			if(table[deltaIn][deltaOut])
				ctrCheck++;
		}
	}


	try{
    GRBModel model = GRBModel(env);

    //State variables
    vector<GRBVar> x(n);
    vector<GRBVar> y(n);
    for(int i = 0; i < n; i++){
		x[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x"+to_string(i));
		y[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y"+to_string(i));
    }
    model.update();

    //SHR constraints
	addSHRRXDiffConstr(model,x,y,beta,gamma);

    model.update();
    model.optimize();

    model.set(GRB_IntParam_PoolSolutions, 1ULL << (2*n));
    model.set(GRB_IntParam_PoolSearchMode, 2);
    model.optimize();

    // Print number of solutions stored
    uint nSolutions = model.get(GRB_IntAttr_SolCount);
    cout << "Number of solutions found: " << nSolutions;
    if(nSolutions == ctrCheck)
    	cout << " : correct number of solutions" << endl;
    else{
    	cout << " : incorrect number of solutions, expected " << ctrCheck << endl;
    	return false;
    }

    for (uint e = 0; e < nSolutions; e++) {
      model.set(GRB_IntParam_SolutionNumber, e);

      vector<uint> valx(n);
      vector<uint> valy(n);
      for(int i = 0; i < n; i++){
      	valx[i] = uint(round(x[i].get(GRB_DoubleAttr_Xn)));
      	valy[i] = uint(round(y[i].get(GRB_DoubleAttr_Xn)));
      }
      uint64_t xint = 0;
      uint64_t yint = 0;
      for(int i = 0; i < n; i++){
      	if(valx[i] != 0) xint |= (1ULL << i);
      	if(valy[i] != 0) yint |= (1ULL << i);
      }
      if(!table[xint][yint]){
      	cout << "Check failed, found solution x = " << xint << " y = " << yint << " with beta = " << beta + " gamma = " << gamma << endl;
      	return false;
      }

    }
    cout << endl;
	return true;

}  catch(GRBException e) {
cout << "Error code = " << e.getErrorCode() << endl;
cout << e.getMessage() << endl;
} catch(...) {
cout << "Exception during optimization" << endl;
}
return false;
}

bool checkCriteria(uint64_t const a,
				   uint64_t const b,
				   uint64_t const d,
				   uint const n,
				   uint const gamma){

	uint64_t mask = (1ULL << n) - 1;
	uint64_t u = (a^b^d) ^ shl(a^b^d,1,n);
	uint64_t v = shl( (a^d) | (b^d) ,1,n);

	// cout << "u   = "; binprint(u,n,true); cout << endl;
	// cout << "v   = "; binprint(v,n,true); cout << endl;
	// cout << "u&v = "; binprint(u&v,n,true); cout << endl;

	if(gamma == 0){
		// cout << "u   = "; binprint(u,n,true); cout << endl;
		// cout << "v   = "; binprint(v,n,true); cout << endl;
		// cout << "u&v = "; binprint(u&v,n,true); cout << endl;
		return ((u&v) == u);
	}
	else{
		uint64_t mask1 = (~(1ULL))&mask;
		uint64_t maskgamma = (~(1ULL << gamma))&mask;
		u = u&mask1&maskgamma;
		v = v&mask1&maskgamma;
		// cout << "u   = "; binprint(u,n,true); cout << endl;
		// cout << "v   = "; binprint(v,n,true); cout << endl;
		// cout << "u&v = "; binprint(u&v,n,true); cout << endl;
		return ((u&v) == u);
	}
}

bool testTheoremModAdd(uint const n,
					   uint const gamma){

	uint64_t bound = (1ULL << n);
	uint64_t mask = bound-1;
	vector<vector<vector<bool>>> ddt(bound, vector<vector<bool>>(bound, vector<bool>(bound,0)));

	for(uint64_t x = 0; x < bound; x++){
		uint64_t rx = CSHL(x,gamma,n);
		for(uint64_t y = 0; y < bound; y++){
			uint64_t z = (x+y)&mask;
			uint64_t ry = CSHL(y,gamma,n);
			uint64_t rz = CSHL(z,gamma,n);
			for(uint64_t dx = 0; dx < bound; dx++){
				uint64_t xp = rx^dx;
				for(uint64_t dy = 0; dy < bound; dy++){
					uint64_t yp = ry^dy;
					uint64_t zp = (xp+yp)&mask;
					uint64_t dz = rz^zp;

					ddt[dx][dy][dz] = true;
				}
			}
		}
	}
	bool check = true;
	for(uint64_t dx = 0; dx < n; dx++){
		for(uint64_t dy = 0; dy < n; dy++){
			for(uint64_t dz = 0; dz < n; dz++){
				if(ddt[dx][dy][dz] != checkCriteria(dx,dy,dz,n,gamma)){
					cout << "n = " << n << " gamma = " << gamma << " dx = " << dx << " dy = " << dy << " dz = " << dz << endl;
					cout << "ddt[dx][dy][dz] = " << ddt[dx][dy][dz] << endl;
					cout << "checkCriteria(dx,dy,dz,n,gamma) = " << checkCriteria(dx,dy,dz,n,gamma) << endl;
					check = false;
				}
			}
		}
	}
	return check;
}



int main(){
	map<bool,uint> ctrTest;

	ctrTest[testXOR2(false)]++;
	ctrTest[testXOR2(true)]++;
	ctrTest[testXOR2(false,1)]++;
	ctrTest[testXOR2(true,1)]++;

	ctrTest[testXOR3(false)]++;
	ctrTest[testXOR3(true)]++;
	ctrTest[testXOR3(false,1)]++;
	ctrTest[testXOR3(true,1)]++;

	ctrTest[testXORn(2)]++;
	ctrTest[testXORn(3)]++;
	ctrTest[testXORn(4)]++;
	ctrTest[testXORn(5)]++;
	ctrTest[testXORn(2,1)]++;
	ctrTest[testXORn(3,1)]++;
	ctrTest[testXORn(4,1)]++;
	ctrTest[testXORn(5,1)]++;

	ctrTest[testRXCstXOR(4)]++;
	ctrTest[testRXCstXOR(7)]++;

	ctrTest[testModAddValue(2)]++;
	ctrTest[testModAddValue(3)]++;
	ctrTest[testModAddValue(4)]++;
	ctrTest[testModAddValue(5)]++;
	ctrTest[testModAddValue(6)]++;

	ctrTest[testSpeckKS(16)]++;
	ctrTest[testSpeckKS(24)]++;
	ctrTest[testSpeckKS(32)]++;
	ctrTest[testSpeckKS(64)]++;

	for(uint n = 3; n < 6; n++){
		for(uint gamma = 0; gamma < n; gamma++)
			ctrTest[testModAddRXDiff(n,gamma)]++;
	}

	ctrTest[testSpeckValue(16,4)]++;
	ctrTest[testSpeckValue(24,3)]++;
	ctrTest[testSpeckValue(64,2)]++;

	ctrTest[testModAddConstValue(2)]++;
	ctrTest[testModAddConstValue(3)]++;
	ctrTest[testModAddConstValue(4)]++;
	ctrTest[testModAddConstValue(5)]++;
	ctrTest[testModAddConstValue(6)]++;

	ctrTest[testLeaKS(128)]++;
	ctrTest[testLeaKS(192)]++;
	ctrTest[testLeaKS(256)]++;

	ctrTest[testLeaValue(128)]++;
	ctrTest[testLeaValue(192)]++;
	ctrTest[testLeaValue(256)]++;

	ctrTest[testRXCstXOR_nonCstVar(3)]++;
	ctrTest[testRXCstXOR_nonCstVar(4)]++;
	ctrTest[testRXCstXOR_nonCstVar(5)]++;

	vector<uint> paramsAlzette({31,24,17,17,0,31,24,16});
	uint32_t cstAlzette = 0xb7e15162;
	ctrTest[testAlzetteValue(paramsAlzette,cstAlzette)]++;
	cstAlzette = 0xbf715880;
	ctrTest[testAlzetteValue(paramsAlzette,cstAlzette)]++;
	cstAlzette = 0x38b4da56;
	ctrTest[testAlzetteValue(paramsAlzette,cstAlzette)]++;
	cstAlzette = 0x324e7738;
	ctrTest[testAlzetteValue(paramsAlzette,cstAlzette)]++;
	cstAlzette = 0xbb1185eb;
	ctrTest[testAlzetteValue(paramsAlzette,cstAlzette)]++;
	cstAlzette = 0x4f7c7b57;
	ctrTest[testAlzetteValue(paramsAlzette,cstAlzette)]++;
	cstAlzette = 0xcfbfa1c8;
	ctrTest[testAlzetteValue(paramsAlzette,cstAlzette)]++;
	cstAlzette = 0xc2b3293d;
	ctrTest[testAlzetteValue(paramsAlzette,cstAlzette)]++;

	vector<uint32_t> keyXTEA({0,0,0,0});
	ctrTest[testXteaKS(keyXTEA)]++;
	keyXTEA = vector<uint32_t>({1,2,3,4});
	ctrTest[testXteaKS(keyXTEA)]++;

	ctrTest[testXteaValue()]++;
	ctrTest[testXteaValue()]++;

	int n = 8;
	for(int alpha = 1; alpha < n; alpha++){
		for(int gamma = 1; gamma < n; gamma++)
			ctrTest[testSHLRXDiff(n,alpha,gamma)]++;
	}

	for(int alpha = 1; alpha < n; alpha++){
		for(int gamma = 1; gamma < n; gamma++)
			ctrTest[testSHRRXDiff(n,alpha,gamma)]++;
	}

	for(uint n = 2; n < 6; n++){
		for(uint gamma = 0; gamma < n; gamma++)
			ctrTest[testTheoremModAdd(n,gamma)]++;
	}


	cout << "------------" << endl;
	cout << ctrTest[true] << " valid tests" << endl;
	cout << ctrTest[false] << " failed tests" << endl;
}
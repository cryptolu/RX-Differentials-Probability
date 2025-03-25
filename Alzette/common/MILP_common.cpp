#include "MILP_common.hpp"
using namespace std;

void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & z,
				  uint const c){
	//Add constraint z = x XOR y XOR c with x,y,z binaries, c a constant (0 or 1)
	if(c == 0){
		m.addConstr(x + y + (1-z) >= 1);
		m.addConstr(x + (1-y) + z >= 1);
		m.addConstr((1-x) + y + z >= 1);
		m.addConstr((1-x) + (1-y) + (1-z) >= 1);
	}
	else{
		m.addConstr(x + y + z >= 1);
		m.addConstr((1-x) + (1-y) + z >= 1);
		m.addConstr((1-x) + y + (1-z) >= 1);
		m.addConstr(x + (1-y) + (1-z) >= 1);
	}
}

void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & z,
				  string const & dumName){
//Add constraint z = x XOR y with x,y,z binaries
//Dummy variable version, probably better not to use
	GRBVar dum = m.addVar(0.0, 1.0, 0.0, GRB_BINARY, dumName);
	m.addConstr(x + y + z == 2*dum);
}

void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & z,
				  uint const c,
				  std::string const & dumName){
//Add constraint z = x XOR y XOR c with x,y,z binaries, c constant
//Dummy variable version, probably better not to use
	if(c == 0)
		addXORConstr(m,x,y,z,dumName);
	else{
		GRBVar dum = m.addVar(0.0, 2.0, 0.0, GRB_INTEGER, dumName);
		m.addConstr(x + y + z + 1 == 2*dum);
	}
}

void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & w,
				  GRBVar & z,
				  uint const c){
	//Add constraint z = x XOR y XOR w XOR c with x,y,w,z binaries, c constant
	if(c == 0){
		m.addConstr(w + x + y + (1-z) >= 1);
		m.addConstr((1-w) + x + y + z >= 1);
		m.addConstr(w + x + (1-y) + z >= 1);
		m.addConstr((1-w) + x + (1-y) + (1-z) >= 1);
		m.addConstr(w + (1-x) + y + z >= 1);
		m.addConstr((1-w) + (1-x) + y + (1-z) >= 1);
		m.addConstr(w + (1-x) + (1-y) + (1-z) >= 1);
		m.addConstr((1-w) + (1-x) + (1-y) + z >= 1);
	}
	else{
		m.addConstr(w + x + y + z >= 1);
		m.addConstr(w + x + (1-y) + (1-z) >= 1);
		m.addConstr(w + (1-x) + y + (1-z) >= 1);
		m.addConstr(w + (1-x) + (1-y) + z >= 1);
		m.addConstr((1-w) + x + y + (1-z) >= 1);
		m.addConstr((1-w) + x + (1-y) + z >= 1);
		m.addConstr((1-w) + (1-x) + y + z >= 1);
		m.addConstr((1-w) + (1-x) + (1-y) + (1-z) >= 1);
	}
}

void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & w,
				  GRBVar & z,
				  std::string const & dumName){
//Add constraint z = x XOR y XOR w with x,y,w,z binaries
//Dummy variable version, probably better not to use
	GRBVar dum = m.addVar(0.0, 2.0, 0.0, GRB_INTEGER, dumName);
	m.addConstr(w + x + y + z == 2*dum);
}

void addXORConstr(GRBModel & m,
				  GRBVar & x,
				  GRBVar & y,
				  GRBVar & w,
				  GRBVar & z,
				  uint const c,
				  std::string const & dumName){
//Add constraint z = x XOR y XOR w XOR c with x,y,w,z binaries, c constant
//Dummy variable version, probably better not to use
	if(c == 0)
		addXORConstr(m,x,y,w,z,c);
	else{
		GRBVar dum = m.addVar(0.0, 2.0, 0.0, GRB_INTEGER, dumName);
		m.addConstr(w + x + y + z + 1 == 2*dum);
	}
}

void addXORConstr(GRBModel & m,
				  std::vector<GRBVar> & vars,
				  GRBVar & z,
				  std::string const & dumName){
//Add constraint z = vars[0] XOR vars[1] XOR ... XOR vars[n-1]
//Arbitrary length version, always use dummy variable if >= 4 input variables
//if 2 or 3 input variables, uses the binary modelization for now as it's likely better

	if(vars.size() == 0)
		cerr << "Error : XOR constraint with 0 input variables, ignored" << endl;
	else if(vars.size() == 1)
		m.addConstr(vars[0] == z);
	else if (vars.size() == 2)
		addXORConstr(m,vars[0], vars[1], z);
	else if (vars.size() == 3)
		addXORConstr(m,vars[0],vars[1],vars[2],z);
	else{
		double bound = ceil((double)(vars.size())/2);
		GRBVar dum = m.addVar(0.0, bound, 0.0, GRB_INTEGER, dumName);
		GRBLinExpr expr = 0;
		for(uint i = 0; i < vars.size(); i++)
			expr += vars[i];
		expr += z;
		m.addConstr(expr == 2*dum);
	}
}
void addXORConstr(GRBModel & m,
				  std::vector<GRBVar> & vars,
				  GRBVar & z,
				  uint const c,
				  std::string const & dumName){
//Add constraint z = vars[0] XOR vars[1] XOR ... XOR vars[n-1] XOR c
//Arbitrary length version, always use dummy variable if >= 4 input variables
//if 2 or 3 input variables, uses the binary modelization for now as it's likely better

	if(c == 0)
		addXORConstr(m,vars,z,dumName);
	else{
		if(vars.size() == 0)
			cerr << "Error : XOR constraint with 0 input variables, ignored" << endl;
		else if(vars.size() == 1)
			m.addConstr(1 - vars[0] == z);
		else if (vars.size() == 2)
			addXORConstr(m,vars[0], vars[1], z,1);
		else if (vars.size() == 3)
			addXORConstr(m,vars[0],vars[1],vars[2],z,1);
		else{
			double bound = ceil((double)(vars.size()+1)/2);
			GRBVar dum = m.addVar(0.0, bound, 0.0, GRB_INTEGER, dumName);
			GRBLinExpr expr = 0;
			for(uint i = 0; i < vars.size(); i++)
				expr += vars[i];
			expr += z;
			m.addConstr(expr + 1 == 2*dum);
		}
	}
}

void addRXCstXORConstr(GRBModel & m,
					   std::vector<GRBVar> & vars,
					   std::vector<GRBVar> & cst,
					   std::vector<GRBVar> & res,
					   uint const s){
//Add a constraint for the RX-diff propagation through res = vars XOR cst where cst is seen as a constant but still a gurobi var (e.g. round key with key-schedule modelized) with RX shift s

	uint n = vars.size();
	if(s > 0){
		for(uint i = 0; i < n; i++)
			addXORConstr(m, vars[i], cst[i], cst[mod(i-s,n)], res[i]);
	}
	else{
		//Regular differential propagation
		for(uint i = 0; i < n; i++)
			m.addConstr(vars[i] == res[i]);
	}
}

void addRXCstXORConstr(GRBModel & m,
					   std::vector<GRBVar> & vars,
					   std::vector<uint> const & cst,
					   std::vector<GRBVar> & res,
					   uint const s){

	uint n = vars.size();
	if(s > 0){
		for(uint i = 0; i < n; i++){
			if((cst[i] ^ cst[mod(i-s,n)]) == 0)
				m.addConstr(vars[i] == res[i]);
			else
				m.addConstr(res[i] == 1 - vars[i]);
		}
	}
	else{
		//Regular differential propagation
		for(uint i = 0; i < n; i++)
			m.addConstr(vars[i] == res[i]);
	}
}

void addSimonCoreConstr(GRBModel & m,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						int a,
						int b,
						int const c,
						std::string const & prefix){
//Add constraints to modelize the differential transition x -> y through the function
//f(x) = (x << a) & (x << b) ^ (x << c)
//x and y are assumed to be of the same size
//Requires the introduction of variables with varname starting by prefix

	int const n = x.size();

	//Requires a > b
	//Rather than checking this, just adjust it if necessary
	if(a <= b){
		uint tmp = a;
		a = b;
		b = tmp;
	}

	
	//Also requires n even and gcd(n,a-b) = 1
	//For safety, make a check
	if(gcd(n,a-b) != 1){
		cerr << "SimonCore constraints with parameters gcd(n,a-b) != 1, not handled" << endl;
		exit(1);
	}
	if(n%2 != 0){
		cerr << "SimonCore constraints with n odd, not handled" << endl;
		exit(1);
	}
	
	vector<GRBVar> gamma(n);
	for(int i = 0; i < n; i++)
		gamma[i] = m.addVar(0.0, 1.0, 0.0, GRB_BINARY, prefix+to_string(i));


	//Indicator constraint to manage x = 11...1 => wt(gamma)%2 == 0
	GRBVar andx = m.addVar(0.0, 1.0, 0.0, GRB_BINARY, prefix+"_andx");
	m.addGenConstrAnd(andx, x.data(), n);
	double bound = ceil((double)(n)/2);
	GRBVar andxdum = m.addVar(0.0, bound, 0.0, GRB_INTEGER, prefix+"_andxdum");
	andxdum.set(GRB_IntAttr_PoolIgnore, 1);
	GRBLinExpr sumgamma = 0;
	for(int i = 0; i < n; i++)
		sumgamma += gamma[i];
	sumgamma -= 2*andxdum;
	m.addGenConstrIndicator(andx,1,sumgamma,GRB_EQUAL,0);

	for(int i = 0; i < n; i++){
		addXORConstr(m,x[mod(i-c,n)],y[i],gamma[i]);
		m.addConstr((1-gamma[i]) + x[mod(i-a,n)] + x[mod(i-b,n)] >= 1);
		m.addConstr(gamma[i] + (1-gamma[mod(i-a+b,n)]) + (1-x[mod(i-b,n)]) + x[mod(i-a,n)] + (1-x[mod(i-2*a+b,n)]) >= 1);
		m.addConstr((1-gamma[i]) + gamma[mod(i-a+b,n)] + (1-x[mod(i-b,n)]) + x[mod(i-a,n)] + (1-x[mod(i-2*a+b,n)]) >= 1);

		// cout << " i = " << i << " (a,b,c) = (" << a << "," << b << "," << c << ")" << endl;
		// cout << "XOR(x" << mod(i-c,n) << ",y" << i << ",gamma" << i << ")" << endl;
		// cout << "gamma" << i << " + (1-x" << mod(i-a,n) << ") + (1-x" << mod(i-b,n) << ") >= 1" << endl;
		// cout << "gamma" << i << " + (1-gamma" << mod(i-a+b,n) << ") + (1-x" << mod(i-b,n) << ") + x" << mod(i-a,n) << " + (1-x" << mod(i-2*a+b,n) << ") >= 1" << endl;
		// cout << "(1-gamma" << i << ") + gamma" << mod(i-a+b,n) << " + (1-x" << mod(i-b,n) << ") + x" << mod(i-a,n) << " + (1-x" << mod(i-2*a+b,n) << ") >= 1" << endl;
		// cout << endl;
	}

}

void addSimonCoreValConstr(GRBModel & model,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						int a,
						int b,
						int const c){
//Add constraints to modelize the transition **in value** y = f(x)
//This version does a direct modelization without intermediate variables

	//f(x) = (x <<< a) & (x <<< b) ^ (x <<< c)
	int const n = x.size();

	for(int i = 0; i < n; i++){
		GRBVar & xa = x[mod(i-a,n)];
		GRBVar & xb = x[mod(i-b,n)];
		GRBVar & xc = x[mod(i-c,n)];
		GRBVar & yi = y[i];
		model.addConstr(xa + xb + xc + (1 - yi) >= 1);
		model.addConstr(xa + xb + (1 - xc) + yi >= 1);
		model.addConstr(xa + (1 - xb) + xc + (1 - yi) >= 1);
		model.addConstr(xa + (1 - xb) + (1 - xc) + yi >= 1);
		model.addConstr((1 - xa) + xb + xc + (1 - yi) >= 1);
		model.addConstr((1 - xa) + xb + (1 - xc) + yi >= 1);
		model.addConstr((1 - xa) + (1 - xb) + xc + yi >= 1);
		model.addConstr((1 - xa) + (1 - xb) + (1 - xc) + (1 - yi) >= 1);
	}
}

void addModAddValueConstr(GRBModel & model, 
						  std::vector<GRBVar> & x,
						  std::vector<GRBVar> & y,
						  std::vector<GRBVar> & z,
						  std::string const & prefix){
//Constraints to represent z = x + y **in value**
//prefix is used to created carry variables

	uint n = x.size();
	//Create carry variables
	vector<GRBVar> c(n-1);
	for(uint i = 0; i < n-1; i++)
		c[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, prefix+"_"+to_string(i));

	//Constraints for the carry
	//c[0] = Maj(x[0], y[0])
	//c[i] = Maj(x[i], y[i], c[i-1]); i > 0

	//i = 0
	//If x[0] or y[0] == 0; then carry[0] = 0
	//Otherwise, carry[0] = 1
	//So this can easily be done with a Min constraint
	GRBVar vars[2] = {x[0], y[0]};
	model.addGenConstrMin(c[0], vars, 2);

	//i > 0
	//Here done with the usual logical modeling
	// x[i]  y[i]  c[i-1]  c[i]
	//  0     0      0      0  Ok
	//  0     0      0      1  Nope
	//  0     0      1      0  Ok
	//  0     0      1      1  Nope
	//  0     1      0      0  Ok
	//  0     1      0      1  Nope
	//  0     1      1      0  Nope
	//  0     1      1      1  Ok
	//  1     0      0      0  Ok
	//  1     0      0      1  Nope
	//  1     0      1      0  Nope
	//  1     0      1      1  Ok
	//  1     1      0      0  Nope
	//  1     1      0      1  Ok
	//  1     1      1      0  Nope
	//  1     1      1      1  Ok
	for(uint i = 1; i < n-1; i++){
		model.addConstr(x[i] + y[i] + c[i-1] + (1 - c[i]) >= 1);             //0  0  0  1 
		model.addConstr(x[i] + y[i] + (1 - c[i-1]) + (1 - c[i]) >= 1);       //0  0  1  1 
		model.addConstr(x[i] + (1 - y[i]) + c[i-1] + (1 - c[i]) >= 1);       //0  1  0  1 
		model.addConstr(x[i] + (1 - y[i]) + (1 - c[i-1]) + c[i] >= 1);       //0  1  1  0 
		model.addConstr((1 - x[i]) + y[i] + c[i-1] + (1 - c[i]) >= 1);       //1  0  0  1 
		model.addConstr((1 - x[i]) + y[i] + (1 - c[i-1]) + c[i] >= 1);       //1  0  1  0 
		model.addConstr((1 - x[i]) + (1 - y[i]) + c[i-1] + c[i] >= 1);       //1  1  0  0 
		model.addConstr((1 - x[i]) + (1 - y[i]) + (1 - c[i-1]) + c[i] >= 1); //1  1  1  0 
	}

	//Constraints for the result
	//z[0] = x[0] ^ y[0]
	//z[i] = x[i] ^ y[i] ^ c[i-1]; i > 0
	addXORConstr(model, x[0], y[0], z[0]);
	for(uint i = 1; i < n; i++)
		addXORConstr(model, x[i], y[i], c[i-1], z[i]);
}

void addModAddConstValueConstr(GRBModel & model, 
							   std::vector<GRBVar> & x,
							   uint64_t const cst,
							   std::vector<GRBVar> & z,
							   std::string const & prefix){
//Constraints to represent z = x + y **in value** where y is a known constant
//prefix is used to create carry variables
//Limited to constants over 64 bits for now

	uint n = x.size();
	//Create carry variables
	vector<GRBVar> c(n-1);
	for(uint i = 0; i < n-1; i++)
		c[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, prefix+"_"+to_string(i));

	//Extract the bits of cst
	vector<uint> y(n,0);
	for(uint i = 0; i < n; i++)
		y[i] = (cst >> i)&1;

	//Constraints for the carry
	//c[0] = Maj(x[0], y[0])
	//c[i] = Maj(x[i], y[i], c[i-1]); i > 0

	//i = 0
	//If x[0] or y[0] == 0; then carry[0] = 0
	//Otherwise, carry[0] = 1
	//So this can easily be done with a Min constraint
	GRBVar vars[1] = {x[0]};
	model.addGenConstrMin(c[0], vars, 1, y[0]);

	//i > 0
	//Done with logical modeling, but we split cases depending on whether y[i] = 0 or 1
	for(uint i = 1; i < n-1; i++){
		if(y[i] == 0){
			model.addConstr(x[i] + c[i-1] + (1 - c[i]) >= 1);             //0  0  0  1 
			model.addConstr(x[i] + (1 - c[i-1]) + (1 - c[i]) >= 1);       //0  0  1  1
			model.addConstr((1 - x[i]) + c[i-1] + (1 - c[i]) >= 1);       //1  0  0  1 
			model.addConstr((1 - x[i]) + (1 - c[i-1]) + c[i] >= 1);       //1  0  1  0 
		}
		else{
			model.addConstr(x[i] + c[i-1] + (1 - c[i]) >= 1);       //0  1  0  1 
			model.addConstr(x[i] + (1 - c[i-1]) + c[i] >= 1);       //0  1  1  0 
			model.addConstr((1 - x[i]) + c[i-1] + c[i] >= 1);       //1  1  0  0 
			model.addConstr((1 - x[i]) + (1 - c[i-1]) + c[i] >= 1); //1  1  1  0 
		}
	}

	//Constraints for the result
	//z[0] = x[0] ^ y[0]
	//z[i] = x[i] ^ y[i] ^ c[i-1]; i > 0
	if(y[0] == 0)
		model.addConstr(z[0] == x[0]);
	else
		model.addConstr(z[0] == (1 - x[0]));

	for(uint i = 1; i < n; i++)
		addXORConstr(model, x[i], c[i-1], z[i], y[i]);
}

void fConstraintRXDiffModAdd(GRBModel & model,
							 GRBVar & x1,
							 GRBVar & x2,
							 GRBVar & x3,
							 GRBVar & x4,
							 GRBVar & x5,
							 GRBVar & x6){
//f-constraint used in the modelization of the RX-diff propagation for mod add

	model.addConstr((1 - x1 ) + x2 + x3 + x4 + x5 + x6 >= 1);
	model.addConstr(x1 + (1 - x2 ) + x3 + x4 + x5 + x6 >= 1);
	model.addConstr(x1 + x2 + (1 - x3 ) + x4 + x5 + x6 >= 1);
	model.addConstr((1 - x1 ) + (1 - x2 ) + (1 - x3 ) + x4 + x5 + x6 >= 1);
	model.addConstr(x1 + (1 - x2) + (1 - x3 ) + (1 - x4 ) + (1 - x5 ) + (1 - x6 ) >= 1);
	model.addConstr((1 - x1 ) + x2 + (1 - x3 ) + (1 - x4 ) + (1 - x5 ) + (1 - x6 ) >= 1);
	model.addConstr((1 - x1 ) + (1 - x2 ) + x3 + (1 - x4 ) + (1 - x5 ) + (1 - x6 ) >= 1);
	model.addConstr(x1 + x2 + x3 + (1 - x4 ) + (1 - x5 ) + (1 - x6 ) >= 1);
}

void addModAddRXDiffConstr(GRBModel & model,
						   std::vector<GRBVar> & a,
						   std::vector<GRBVar> & b,
						   std::vector<GRBVar> & d,
						   uint const gamma){
//Add constraints to modelize the RX-diff propagation for mod add
// (a,gamma),(b,gamma) -> (d,gamma)

	uint const n = a.size();

	//constraint on the first bit only if gamma == 0
	if(gamma == 0)
		addXORConstr(model, a[0], b[0], d[0]);

	for(uint i = 1; i < gamma; i++) //noop if gamma==0, as expected
		fConstraintRXDiffModAdd(model,a[i],b[i],d[i],a[i-1],b[i-1],d[i-1]);
		
	//no constraint for i = gamma if gamma != 0
	//if gamma==0, the next loop adds all the cases for 1 <= i < n

	for(uint i = gamma+1; i < n; i++)
		fConstraintRXDiffModAdd(model,a[i],b[i],d[i],a[i-1],b[i-1],d[i-1]);

}


void addSHLRXDiffConstr(GRBModel & model,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						int const alpha,
						int const gamma){
//Add constraints to modelize the RX-Diff propagation for y = x << alpha

	int n = x.size();

	if(gamma+alpha <= n){
		if(gamma <= alpha){
			//0 <= i < gamma --> anything is valid
			//No constraints

			//gamma <= i < alpha --> y[i] = 0
			for(int i = gamma; i < alpha; i++)
				model.addConstr(y[i] == 0);

			//alpha <= i < gamma+alpha --> anything is valid
			//No constraints
			
			//gamma+alpha <= i < n --> y[i] = x[i-alpha]
			for(int i = gamma+alpha; i < n; i++)
				model.addConstr(y[i] == x[mod(i-alpha,n)]);
		}
		else{ //gamma > alpha
			//0 <= i < alpha --> anything is valid
			//No constraints

			//alpha <= i < gamma --> y[i] = x[i-alpha]
			for(int i = alpha; i < gamma; i++)
				model.addConstr(y[i] == x[mod(i-alpha,n)]);

			//gamma <= i < gamma+alpha --> anything is valid
			//No constraints

			//gamma+alpha <= i < n --> y[i] = x[i-alpha]
			for(int i = gamma+alpha; i < n; i++)
				model.addConstr(y[i] == x[mod(i-alpha,n)]);
		}
	}
	else{ //gamma+alpha > n
		if(gamma <= alpha){
			//0 <= i < gamma+alpha-n --> y[i] = 0
			for(int i = 0; i < gamma+alpha-n; i++)
				model.addConstr(y[i] == 0);

			//gamma+alpha-n <= i < gamma --> anything is valid
			//No constraints

			//gamma <= i < alpha --> y[i] = 0
			for(int i = gamma; i < alpha; i++)
				model.addConstr(y[i] == 0);

			//alpha <= i < n --> anything is valid
			//No constraints
		}
		else{ //gamma > alpha
			//0 <= i < gamma+alpha-n --> y[i] = 0
			for(int i = 0; i < gamma+alpha-n; i++)
				model.addConstr(y[i] == 0);

			//gamma+alpha-n <= i < alpha --> anything is valid
			//No constraints

			//alpha <= i < gamma --> y[i] = x[i-alpha]
			for(int i = alpha; i < gamma; i++)
				model.addConstr(y[i] == x[mod(i-alpha,n)]);

			//gamma <= i < n --> anything is valid
			//No constraints
		}
	}
}

void addSHRRXDiffConstr(GRBModel & model,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						int const beta,
						int const gamma){
//Add constraints to modelize the RX-Diff propagation for y = x << beta

	int n = x.size();
	
	if(gamma+beta <= n){
		if(gamma <= beta){
			//0 <= i < gamma --> anything is valid
			//No constraints

			//gamma <= i < n-beta --> y[i] = x[i+beta]
			for(int i = gamma; i < n-beta; i++)
				model.addConstr(y[i] == x[mod(i+beta,n)]);

			//n-beta <= i < n-beta+gamma --> anything is valid
			//No constraints

			//n-beta+gamma <= i < n --> y[i] = 0
			for(int i = n-beta+gamma; i < n; i++)
				model.addConstr(y[i] == 0);
		}
		else{ //gamma > beta
			//0 <= i < gamma-beta --> y[i] = x[i+beta]
			for(int i = 0; i < gamma-beta; i++)
				model.addConstr(y[i] == x[mod(i+beta,n)]);

			//gamma-beta <= i < gamma --> anything is valid
			//No constraints

			//gamma <= i < n-beta --> y[i] = x[i+beta]
			for(int i = gamma; i < n-beta; i++)
				model.addConstr(y[i] == x[mod(i+beta,n)]);

			//n-beta <= i < n --> anything is valid
			//No constraints
		}
	}
	else{ //gamma+beta > n
		if(gamma <= beta){
			//0 <= i < n-beta --> anything is valid
			//No constraints

			//n-beta <= i < gamma --> y[i] = 0
			for(int i = n-beta; i < gamma; i++)
				model.addConstr(y[i] == 0);

			//gamma <= i < gamma+n-beta --> anything is valid
			//No constraints

			//gamma+n-beta <= i < n --> y[i] = 0
			for(int i = gamma+n-beta; i < n; i++)
				model.addConstr(y[i] == 0);
		}
		else{ //gamma > beta
			//0 <= i < gamma-beta --> y[i] = x[i+beta]
			for(int i = 0; i < gamma-beta; i++)
				model.addConstr(y[i] == x[mod(i+beta,n)]);

			//gamma-beta <= i < n-beta --> anything is valid
			//No constraints

			//n-beta <= i < gamma --> y[i] = 0
			for(int i = n-beta; i < gamma; i++)
				model.addConstr(y[i] == 0);

			//gamma <= i < n --> anything is valid
			//No constraints
		}
	}
}

void addAllEqualBoolConstr(GRBModel & model,
						   std::vector<GRBVar> & x,
						   GRBVar & eq){
//Add a constraints so that eq == 1 iif all variables in x have the same value

	int n = x.size();
	// # eq = 1 =>  xs = all equal  (via x1 + x2 + ... + x_(n-1) - (n-1)*x = 0)
    GRBLinExpr sumx = 0;
    for(int i = 0; i < n-1; i++)
    	sumx += x[i];
    model.addConstr( sumx - (n-1)*x[n-1] + (n-1)*(1 - eq) >= 0);
    model.addConstr(-sumx + (n-1)*x[n-1] + (n-1)*(1 - eq) >= 0);

    // # eq = 0 => non-eq  (via 1 <= sum(xs) <= n-1)
    model.addConstr( sumx + x[n-1] + eq >= 1);
    model.addConstr(-sumx - x[n-1] + eq >= -n+1);
}

void addNotAllEqualBoolConstr(GRBModel & model,
							  std::vector<GRBVar> & x,
							  GRBVar & eq){
//Add a constraints so that eq == 0 iif all variables in x have the same value
	//Essentially, same as above but replacing eq by 1-eq
	int n = x.size();
	// # eq = 1 =>  xs = all equal  (via x1 + x2 + ... + x_(n-1) - (n-1)*x = 0)
    GRBLinExpr sumx = 0;
    for(int i = 0; i < n-1; i++)
    	sumx += x[i];
    model.addConstr( sumx - (n-1)*x[n-1] + (n-1)*eq >= 0);
    model.addConstr(-sumx + (n-1)*x[n-1] + (n-1)*eq >= 0);

    // # eq = 0 => non-eq  (via 1 <= sum(xs) <= n-1)
    model.addConstr( sumx + x[n-1] + (1-eq) >= 1);
    model.addConstr(-sumx - x[n-1] + (1-eq) >= -n+1);
}

void modelMux(GRBModel & model,
			  GRBVar & dst,
			  GRBVar & flag,
			  GRBLinExpr & expr0,
			  GRBLinExpr & expr1,
			  double const M){
	//Combined together, these constraints modelize
	//If flag == 0, then dst == expr0
	//If flag == 1, then dst == expr1
	model.addConstr(dst >= expr0 - M * flag); //If flag == 0, then dst >= expr0
    model.addConstr(dst <= expr0 + M * flag); //If flag == 0, then dst <= expr0
    model.addConstr(dst >= expr1 - M * (1-flag)); //If flag == 1, then dst >= expr1
    model.addConstr(dst <= expr1 + M * (1-flag)); //If flag == 1, then dst <= expr1
}

void modelMux(GRBModel & model,
			  GRBVar & dst,
			  GRBVar & flag,
			  GRBLinExpr & expr0,
			  GRBLinExpr & expr1){
    //Or using Indicator constraints (not sure which one is better, need to test)
    model.addGenConstrIndicator(flag,0,dst-expr0,GRB_EQUAL,0);
    model.addGenConstrIndicator(flag,1,dst-expr1,GRB_EQUAL,0);
}

void addTnConstr(GRBModel & model,
				 std::vector<GRBVar> & alpha,
				 std::vector<GRBVar> & beta,
				 std::vector<GRBVar> & delta,
				 GRBVar & lsbOther,
				 GRBVar & wt,
				 std::string const & suffix){
	/*
	Add constraints to the model to modelize the relation
	wt = -log2(Tn(alpha,beta,delta,lsbOther))
	suffix is used to create intermediary variables with unique names
	n is obtained from the size of the alpha vector
	*/

	int n = alpha.size();

	if(n > 1){
		//Create intermediary variables
		vector<GRBVar> xorvars(n);
		for(int i = 0; i < n; i++)
			xorvars[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "xorvars"+to_string(i)+suffix);

		vector<GRBVar> neqvars(n-1);
		for(int i = 0; i < n-1; i++)
			neqvars[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "neqvars"+to_string(i)+suffix);
		
		GRBVar d = model.addVar(0,n-1,0,GRB_INTEGER, "d"+suffix);
		GRBVar xorAllSame = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "xorAllSame"+suffix);


		//Compute the log tables
		//Technically, could be precomputed, but the cost is very negligeable
		vector<double> table_log_pos(n);
		vector<double> table_log_neg(n);
		for(int i = 0; i < n; i++){
			table_log_pos[i] = log2(1 + pow(2,-n+i));
			table_log_neg[i] = log2(1 - pow(2,-n+i));
		}

		auto const [p_minpos,p_maxpos] = minmax_element(begin(table_log_pos), end(table_log_pos));
		auto const [p_minneg,p_maxneg] = minmax_element(begin(table_log_neg), end(table_log_neg));
		auto minpos = *p_minpos;
		auto maxpos = *p_maxpos;
		auto minneg = *p_minneg;
		auto maxneg = *p_maxneg;

		//log variables
		double eps = 0.001;
		GRBVar lg_pos = model.addVar(minpos-eps, maxpos+eps,0,GRB_CONTINUOUS,"lgPos"+suffix);
		GRBVar lg_neg = model.addVar(minneg-eps, maxneg+eps,0,GRB_CONTINUOUS,"lgNeg"+suffix);
		GRBVar lg = model.addVar(minneg-eps, maxpos+eps,0,GRB_CONTINUOUS,"lg"+suffix);

		//Xor constraints
		for(int i = 0; i < n; i++)
			addXORConstr(model, alpha[i], beta[i], delta[i], xorvars[i]);

		//neq and d
		for(int i = 0; i < n-1; i++){
			vector<GRBVar> vars({alpha[i],beta[i],delta[i]});
			addNotAllEqualBoolConstr(model, vars, neqvars[i]);
		}
		GRBLinExpr sum_neq = 0;
		for(int i = 0; i < n-1; i++)
			sum_neq += neqvars[i];
		model.addConstr(sum_neq == d);

		//xorAllSame
		addAllEqualBoolConstr(model,xorvars,xorAllSame);

		//lg_pos, lg_neg
		vector<double> xpts(n);
		for(int i = 0; i < n; i++)
			xpts[i] = i;
		model.addGenConstrPWL(d, lg_pos, n, xpts.data(), table_log_pos.data());
		model.addGenConstrPWL(d, lg_neg, n, xpts.data(), table_log_neg.data());

		//lg
		GRBLinExpr expr_lg_pos = lg_pos;
		GRBLinExpr expr_lg_neg = lg_neg;
		modelMux(model,lg,lsbOther,expr_lg_pos,expr_lg_neg);

		//wt
		GRBLinExpr expr_dp1 = d+1;
		GRBLinExpr expr_dp1mlg = d+1-lg;
		modelMux(model,wt,xorAllSame,expr_dp1,expr_dp1mlg);

		model.update();
	}
	else{ //n = 1
		//This case is simpler
		//Compared to above, we always have xorAllSame == 1 and d == 0
		//We also don't need to use any PWL, the probability only has 2 values depending on lsbOther
		double lg_pos = 1 - log2(1 + pow(2,-n));
		double lg_neg = 1 - log2(1 - pow(2,-n));
		//lsbOther == 0 -> wt = lg_pos
		//lsbOther == 1 -> wt = lg_neg
		GRBLinExpr expr_lg_pos = lg_pos;
		GRBLinExpr expr_lg_neg = lg_neg;
		modelMux(model,wt,lsbOther,expr_lg_pos,expr_lg_neg);

		model.update();

	}
}

void addModAddRXProbaConstr(GRBModel & model,
							std::vector<GRBVar> & alpha,
							std::vector<GRBVar> & beta,
							std::vector<GRBVar> & delta,
							GRBVar & wt,
							int const k,
							std::string const & suffix){
	/*
	Add constraints to modelize the probability of the RX-differential (alpha,beta) -> delta
	through the modular addition, with RX-rotation k
	The result is given by wt = -log2(probability)
	suffix is used to create intermediary variables with unique names
	n is obtained from the size of the alpha vector
	*/

	int n = alpha.size();

	//Make 2 variables for the intermediate weights for each half
	GRBVar wtL = model.addVar(0,2*(n-k),0,GRB_CONTINUOUS, "wtL"+suffix);
	GRBVar wtR = model.addVar(0,2*k,0,GRB_CONTINUOUS, "wtR"+suffix);

	//Make 2 variables the XOR of the LSBs of each half
	GRBVar lsbL = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "lsbL"+suffix);
	GRBVar lsbR = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "lsbR"+suffix);

	model.update();

	//Add the constraints for the XOR of the LSBs
	addXORConstr(model, alpha[0], beta[0], delta[0], lsbR);
	addXORConstr(model, alpha[k], beta[k], delta[k], lsbL);

	//Tn constraints for the right half
	vector<GRBVar> alphaR(k);
	vector<GRBVar> betaR(k);
	vector<GRBVar> deltaR(k);
	for(int i = 0; i < k; i++){
		alphaR[i] = alpha[i];
		betaR[i] = beta[i];
		deltaR[i] = delta[i];
	}
	addTnConstr(model, alphaR, betaR, deltaR, lsbL, wtR, "R"+suffix);

	//Tn constraints for the left half
	vector<GRBVar> alphaL(n-k);
	vector<GRBVar> betaL(n-k);
	vector<GRBVar> deltaL(n-k);
	for(int i = 0; i < (n-k); i++){
		alphaL[i] = alpha[k+i];
		betaL[i] = beta[k+i];
		deltaL[i] = delta[k+i];
	}
	addTnConstr(model, alphaL, betaL, deltaL, lsbR, wtL, "L"+suffix);

	//Final weight is the sum of the two intermediate weights
	model.addConstr(wt == wtL + wtR);
}

void addVectorEqualValConstr(GRBModel & model,
							 std::vector<GRBVar> & vars,
							 std::vector<uint> const & val){
/*
Add constraints so that vars[i] == val[i]
*/
	for(uint i = 0; i < vars.size(); i++){
		model.addConstr(vars[i] == val[i]);
	}
}



uint64_t getUintSolutionFromBinVector(std::vector<GRBVar> & vars){
//Given a vector of binary variables as input, return the corresponding integer
	uint64_t res = 0;
	for(uint i = 0; i < vars.size(); i++){
		uint64_t val = uint64_t(round(vars[i].get(GRB_DoubleAttr_X)));
		res |= (val << i);
	}
	return res;
}


std::vector<GRBVar> getArrayVar(GRBModel & model,
								uint const dim,
								std::string const & prefix){
//Get all vars in an array s.t. t[i] is the variable names prefix+i
	vector<GRBVar> res(dim);
	for(uint i = 0; i < dim; i++)
		res[i] = model.getVarByName(prefix+to_string(i));
	return res;
}

std::vector<std::vector<GRBVar>> 
get2DArrayVar(GRBModel & model,
			  uint const dim1,
			  uint const dim2,
			  std::string const & prefix){
//Get all vars in a 2D array s.t. t[i][j] is the variable names prefix+i+_+j
	vector<vector<GRBVar>> res(dim1, vector<GRBVar>(dim2));
	for(uint i = 0; i < dim1; i++){
		res[i] = getArrayVar(model, dim2, prefix+to_string(i)+"_");
	}
	return res;
}

std::vector<GRBVar> genArrayBinVar(GRBModel & model,
								   uint const dim,
								   std::string const & prefix){
//Return a vector of binary variables named prefix+i
	vector<GRBVar> res(dim);
	for(uint i = 0; i < dim; i++)
		res[i] = model.addVar(0.0,1.0,0.0,GRB_BINARY,prefix+to_string(i));
	return res;
}	

std::vector<std::vector<GRBVar>> genArray2DBinVar(GRBModel & model,
												  uint const dim1,
												  uint const dim2,
												  std::string const & prefix){
//Return a 2D vector of binary variables named prefix+i+_+j
	vector<vector<GRBVar>> res(dim1);
	for(uint i = 0; i < dim1; i++)
		res[i] = genArrayBinVar(model,dim2,prefix+to_string(i)+"_");
	return res;
}

void addRXDifferentialValueBinding(GRBModel & model,
								   std::vector<GRBVar> & x0,
								   std::vector<GRBVar> & x1,
								   std::vector<GRBVar> & dx,
								   uint const gamma){
//Add the constraints to enforce dx = rotl(x0,gamma) xor x1
	uint n = x0.size();
	vector<GRBVar> tmp(n);
	for(uint i = 0; i < n; i++)
		tmp[i] = x0[mod(i-gamma,n)]; //x0 rotated to the left by gamma
	for(uint i = 0; i < n; i++)
		addXORConstr(model, tmp[i], x1[i], dx[i]);
}

void addValueToBinVectorConstr(GRBModel & model,
							   std::vector<GRBVar> & vars,
							   uint64_t val){
//add constraints vars[i] == val[i] (val[i] = i-th bit of val)
	uint n = vars.size();
	for(uint i = 0; i < n; i++){
		model.addConstr(vars[i] == (val&1));
		val >>= 1;
	}
}

void
addMatsuiLikeConstr(GRBModel & model,
					std::map<uint,uint> const & lowerBoundConsecutiveRounds,
					std::vector<GRBVar> & vars,
					GRBVar & sumneq){
//Add matsui-like constraints to derive lower bounds on consecutive rounds
//e.g. if lowerBoundConsecutiveRounds[3] = 5, then any sum of three consecutive variables in vars must be >= 5
//Also add constraints from ZSCH2018 but those might be redundant

	uint nbRound = vars.size();
    cout << "Currently known bounds on less rounds:" << endl;
    for(auto const & round_bound : lowerBoundConsecutiveRounds){
    	auto const rnd = round_bound.first;
    	auto const bound = round_bound.second;
    	cout << rnd << " rounds, >= " << bound << endl;

    	//Basic bounds on consecutive rounds
    	cout << "Constraints on consecutive rounds:";
    	for(uint start = 0; start <= (nbRound-rnd); start++){
			GRBLinExpr sumconsec = 0;
			cout << "sum over rounds ";
			for(uint r = start; r < start+rnd; r++){
				cout << r << " ";
				sumconsec += vars[r];
			}
			cout << " >= " << bound << endl;
			model.addConstr(sumconsec >= bound);
    	}

    	//Matsui-like bounds from ZSCH2018 (redundant with bounds above ?)
    	//sum(first (nbRounds-rnd) rounds) + Best(rnd rounds) <= sumneq
    	cout << "Matsui-like 1: sum over rounds ";
    	GRBLinExpr matsuiLike1 = 0;
    	for(uint r = 0; r < (nbRound-rnd); r++){
    		cout << r << " ";
    		matsuiLike1 += vars[r];
    	}
    	cout << " + " << bound << " <= sumneq" << endl;
    	model.addConstr(matsuiLike1 + bound <= sumneq);

    	//sum(last (nbRounds-rnd) rounds) + Best(rnd rounds) <= sumneq
    	GRBLinExpr matsuiLike2 = 0;
    	cout << "Matsui-like 2: sum over rounds ";
    	for(uint r = rnd; r < nbRound; r++){
    		cout << r << " ";
    		matsuiLike2 += vars[r];
    	}
    	cout << " + " << bound << " <= sumneq" << endl;
    	model.addConstr(matsuiLike2 + bound <= sumneq);
    }

}
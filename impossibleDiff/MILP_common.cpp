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
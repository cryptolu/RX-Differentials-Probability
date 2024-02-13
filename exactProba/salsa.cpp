#include "salsa.hpp"
using namespace std;

std::vector<GRBVar>
addSalsaOpConstr(GRBModel & model,
				 vector<GRBVar> & a,
				 vector<GRBVar> & b,
				 vector<GRBVar> & c,
				 vector<GRBVar> & out,
				 uint const rot,
				 uint const gamma,
				 std::string const & prefix){
/*
Add the constraints to modelize RX-differentials of out = ((a + b) <<< rot) ^ c
Create temp variables named prefix+"_"+i for i in [0,31]
Return the vector of temp vars
*/

	//First, tmp = (a + b)
	auto tmp = genArrayBinVar(model,32,prefix);
	addModAddRXDiffConstr(model,a,b,tmp,gamma);
	//Rotate tmp
	vector<GRBVar> rotTmp(32);
	for(uint i = 0; i < 32; i++)
		rotTmp[i] = tmp[mod(i-rot,32)];
	//out = rotTmp ^ c
	for(uint i = 0; i < 32; i++)
		addXORConstr(model,rotTmp[i],c[i],out[i]);
	return tmp;
}

std::vector<std::vector<GRBVar>>
addSalsaQuarterRoundConstr(GRBModel & model,
						   uint const gamma,
						   std::array<std::vector<GRBVar>,4> & in,
						   std::array<std::vector<GRBVar>,4> & out,
						   std::string const & suffix){
/*
Add constraints between in and out to modelize the RXDiff propagation from in to out
in and out should be of size (4,32)
suffix is used for unique variable names
Return an array of 4 vectors of variables used to modelize the mod add
Quarter round is defined as 
input {a,b,c,d}
b ^= (a + d) <<< 7;
c ^= (b + a) <<< 9;
d ^= (c + b) <<< 13;
a ^= (d + c) <<< 18;
output {a,b,c,d}
*/
	std::vector<std::vector<GRBVar>> z(4);
	//b ^= (a + d) <<< 7 i.e. out[1] = ((in[0] + in[3]) <<< 7) ^ in[1]
	z[0] = addSalsaOpConstr(model,in[0],in[3],in[1],out[1],7,gamma,"dz"+suffix+"_0_");

	//c ^= (b + a) <<< 9 i.e. out[2] = ((out[1] + in[0]) <<< 9) ^ in[2]
	z[1] = addSalsaOpConstr(model,out[1],in[0],in[2],out[2],9,gamma,"dz"+suffix+"_1_");

	//d ^= (c + b) <<< 13 i.e. out[3] = ((out[2] + out[1]) <<< 13) ^ in[3]
	z[2] = addSalsaOpConstr(model,out[2],out[1],in[3],out[3],13,gamma,"dz"+suffix+"_2_");

	//a ^= (d + c) <<< 18 i.e. out[0] = ((out[3] + out[2]) <<< 18) ^ in[0]
	z[3] = addSalsaOpConstr(model,out[3],out[2],in[0],out[0],18,gamma,"dz"+suffix+"_3_");
	return z;
}

std::pair<std::vector<std::vector<std::vector<GRBVar>>>, 
		  std::vector<std::vector<std::vector<std::vector<GRBVar>>>>>
addSalsaRXDConstr(GRBModel & model,
				  uint const nbRound,
				  uint const startingRound,
				  uint const gamma){
/*
Add constraints for RX-differential propagation for Salsa
Return the variables {x,z} used in the modelization
*/

	uint n = 32;
	vector<vector<vector<GRBVar>>> x(nbRound+1);
	vector<vector<vector<vector<GRBVar>>>> z(nbRound, vector<vector<vector<GRBVar>>>(4));
	for(uint r = 0; r < nbRound+1; r++)
		x[r] = genArray2DBinVar(model, 16, n, "dx"+to_string(r)+"_");

	for(uint r = startingRound; r < startingRound+nbRound; r++){
		uint i = r-startingRound;
		auto &  xr = x[i];
		auto & xr1 = x[i+1];
		std::array<std::vector<GRBVar>,4>  in;
		std::array<std::vector<GRBVar>,4> out;

		if(r%2 == 0){
			//QR(x[ 0], x[ 4], x[ 8], x[12]);	// column 1
			 in[0] =  xr[0];  in[1] =  xr[4];  in[2] =  xr[8];  in[3] =  xr[12];
			out[0] = xr1[0]; out[1] = xr1[4]; out[2] = xr1[8]; out[3] = xr1[12];
			z[i][0] = addSalsaQuarterRoundConstr(model,gamma,in,out, to_string(i)+"_0");

			//QR(x[ 5], x[ 9], x[13], x[ 1]);	// column 2
			 in[0] =  xr[5];  in[1] =  xr[9];  in[2] =  xr[13];  in[3] =  xr[1];
			out[0] = xr1[5]; out[1] = xr1[9]; out[2] = xr1[13]; out[3] = xr1[1];
			z[i][1] = addSalsaQuarterRoundConstr(model,gamma,in,out, to_string(i)+"_1");

			//QR(x[10], x[14], x[ 2], x[ 6]);	// column 3
			 in[0] =  xr[10];  in[1] =  xr[14];  in[2] =  xr[2];  in[3] =  xr[6];
			out[0] = xr1[10]; out[1] = xr1[14]; out[2] = xr1[2]; out[3] = xr1[6];
			z[i][2] = addSalsaQuarterRoundConstr(model,gamma,in,out, to_string(i)+"_2");

			//QR(x[15], x[ 3], x[ 7], x[11]);	// column 4
			 in[0] =  xr[15];  in[1] =  xr[3];  in[2] =  xr[7];  in[3] =  xr[11];
			out[0] = xr1[15]; out[1] = xr1[3]; out[2] = xr1[7]; out[3] = xr1[11];
			z[i][3] = addSalsaQuarterRoundConstr(model,gamma,in,out, to_string(i)+"_3");
		}
		else{
			//QR(x[ 0], x[ 1], x[ 2], x[ 3]);	// row 1
			 in[0] =  xr[0];  in[1] =  xr[1];  in[2] =  xr[2];  in[3] =  xr[3];
			out[0] = xr1[0]; out[1] = xr1[1]; out[2] = xr1[2]; out[3] = xr1[3];
			z[i][0] = addSalsaQuarterRoundConstr(model,gamma,in,out, to_string(i)+"_0");

			//QR(x[ 5], x[ 6], x[ 7], x[ 4]);	// row 2
			 in[0] =  xr[5];  in[1] =  xr[6];  in[2] =  xr[7];  in[3] =  xr[4];
			out[0] = xr1[5]; out[1] = xr1[6]; out[2] = xr1[7]; out[3] = xr1[4];
			z[i][1] = addSalsaQuarterRoundConstr(model,gamma,in,out, to_string(i)+"_1");

			//QR(x[10], x[11], x[ 8], x[ 9]);	// row 3
			 in[0] =  xr[10];  in[1] =  xr[11];  in[2] =  xr[8];  in[3] =  xr[9];
			out[0] = xr1[10]; out[1] = xr1[11]; out[2] = xr1[8]; out[3] = xr1[9];
			z[i][2] = addSalsaQuarterRoundConstr(model,gamma,in,out, to_string(i)+"_2");

			//QR(x[15], x[12], x[13], x[14]);	// row 4
			 in[0] =  xr[15];  in[1] =  xr[12];  in[2] =  xr[13];  in[3] =  xr[14];
			out[0] = xr1[15]; out[1] = xr1[12]; out[2] = xr1[13]; out[3] = xr1[14];
			z[i][3] = addSalsaQuarterRoundConstr(model,gamma,in,out, to_string(i)+"_3");
		}
	}
	return make_pair(x,z);
}

void addNeqConstrSingle(GRBModel & model,
						std::vector<GRBVar> & x,
						std::vector<GRBVar> & y,
						std::vector<GRBVar> & z,
						std::vector<GRBVar> & neq,
						uint const gamma){
	//Add constraints neq[i] = NotAllEqual(x[i],y[i],z[i]) for all i != (gamma-1)
	for(uint i = 0; i < gamma-1; i++){
		vector<GRBVar> vars({x[i],y[i],z[i]});
		addNotAllEqualBoolConstr(model,vars,neq[i]);
	}
	for(uint i = gamma; i < neq.size(); i++){
		vector<GRBVar> vars({x[i],y[i],z[i]});
		addNotAllEqualBoolConstr(model,vars,neq[i]);
	}
}

void 
addSalsaQuarterRoundNeqConstr(GRBModel & model,
							  std::array<std::vector<GRBVar>,4> &  in,
							  std::array<std::vector<GRBVar>,4> & out,
							  std::vector<std::vector<GRBVar>> & z,
							  std::vector<std::vector<GRBVar>> & neq,
							  uint const gamma){

	addNeqConstrSingle(model, in[0], in[3],z[0], neq[0],gamma);
	addNeqConstrSingle(model,out[1], in[0],z[1], neq[1],gamma);
	addNeqConstrSingle(model,out[2],out[1],z[2], neq[2],gamma);
	addNeqConstrSingle(model,out[3],out[2],z[3], neq[3],gamma);
}

void 
addSalsaNeqConstr(GRBModel & model,
				  uint const nbRound,
				  uint const startingRound,
				  uint const gamma,
				  std::vector<std::vector<std::vector<GRBVar>>> & x,
				  std::vector<std::vector<std::vector<std::vector<GRBVar>>>> & z,
				  std::vector<std::vector<std::vector<std::vector<GRBVar>>>> & neq){
	//Link neq vars to state vars

	for(uint r = startingRound; r < startingRound+nbRound; r++){
		uint i = r-startingRound;
		auto &  xr = x[i];
		auto & xr1 = x[i+1];
		std::array<std::vector<GRBVar>,4>  in;
		std::array<std::vector<GRBVar>,4> out;

		if(r%2 == 0){
			//QR(x[ 0], x[ 4], x[ 8], x[12]);	// column 1
			 in[0] =  xr[0];  in[1] =  xr[4];  in[2] =  xr[8];  in[3] =  xr[12];
			out[0] = xr1[0]; out[1] = xr1[4]; out[2] = xr1[8]; out[3] = xr1[12];
			addSalsaQuarterRoundNeqConstr(model,in,out,z[i][0],neq[i][0],gamma);	

			//QR(x[ 5], x[ 9], x[13], x[ 1]);	// column 2
			 in[0] =  xr[5];  in[1] =  xr[9];  in[2] =  xr[13];  in[3] =  xr[1];
			out[0] = xr1[5]; out[1] = xr1[9]; out[2] = xr1[13]; out[3] = xr1[1];
			addSalsaQuarterRoundNeqConstr(model,in,out,z[i][1],neq[i][1],gamma);

			//QR(x[10], x[14], x[ 2], x[ 6]);	// column 3
			 in[0] =  xr[10];  in[1] =  xr[14];  in[2] =  xr[2];  in[3] =  xr[6];
			out[0] = xr1[10]; out[1] = xr1[14]; out[2] = xr1[2]; out[3] = xr1[6];
			addSalsaQuarterRoundNeqConstr(model,in,out,z[i][2],neq[i][2],gamma);

			//QR(x[15], x[ 3], x[ 7], x[11]);	// column 4
			 in[0] =  xr[15];  in[1] =  xr[3];  in[2] =  xr[7];  in[3] =  xr[11];
			out[0] = xr1[15]; out[1] = xr1[3]; out[2] = xr1[7]; out[3] = xr1[11];
			addSalsaQuarterRoundNeqConstr(model,in,out,z[i][3],neq[i][3],gamma);
		}
		else{
			//QR(x[ 0], x[ 1], x[ 2], x[ 3]);	// row 1
			 in[0] =  xr[0];  in[1] =  xr[1];  in[2] =  xr[2];  in[3] =  xr[3];
			out[0] = xr1[0]; out[1] = xr1[1]; out[2] = xr1[2]; out[3] = xr1[3];
			addSalsaQuarterRoundNeqConstr(model,in,out,z[i][0],neq[i][0],gamma);

			//QR(x[ 5], x[ 6], x[ 7], x[ 4]);	// row 2
			 in[0] =  xr[5];  in[1] =  xr[6];  in[2] =  xr[7];  in[3] =  xr[4];
			out[0] = xr1[5]; out[1] = xr1[6]; out[2] = xr1[7]; out[3] = xr1[4];
			addSalsaQuarterRoundNeqConstr(model,in,out,z[i][1],neq[i][1],gamma);

			//QR(x[10], x[11], x[ 8], x[ 9]);	// row 3
			 in[0] =  xr[10];  in[1] =  xr[11];  in[2] =  xr[8];  in[3] =  xr[9];
			out[0] = xr1[10]; out[1] = xr1[11]; out[2] = xr1[8]; out[3] = xr1[9];
			addSalsaQuarterRoundNeqConstr(model,in,out,z[i][2],neq[i][2],gamma);

			//QR(x[15], x[12], x[13], x[14]);	// row 4
			 in[0] =  xr[15];  in[1] =  xr[12];  in[2] =  xr[13];  in[3] =  xr[14];
			out[0] = xr1[15]; out[1] = xr1[12]; out[2] = xr1[13]; out[3] = xr1[14];
			addSalsaQuarterRoundNeqConstr(model,in,out,z[i][3],neq[i][3],gamma);
		}
	}
}

GRBModel
getSalsaModel(uint const nbRound,
			  uint const startingRound,
			  uint const gamma,
			  GRBEnv & env){
/*
Generate a model of nbRound rounds of Salsa, starting from startingRound
startingRound can be used to adjust whether the first round should be even or odd
*/

	uint n = 32;
	GRBModel model = GRBModel(env);
	auto [x,z] = addSalsaRXDConstr(model,nbRound,startingRound,gamma);

	//Neq vars
	vector<vector<vector<vector<GRBVar>>>> neq(nbRound,
		   vector<vector<vector<GRBVar>>>	  (4,
				  vector<vector<GRBVar>>	  (4,
						 vector<GRBVar>		  (n-1))));

	for(uint r = 0; r < nbRound; r++){
		for(uint i = 0; i < 4; i++){
			for(uint j = 0; j < 4; j++){
				for(uint k = 0; k < gamma-1; k++)
					neq[r][i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, 
						"neq"+to_string(r)+"_"+to_string(i)+"_"+to_string(j)+"_"+to_string(k));
				for(uint k = gamma; k < n-1; k++)
					neq[r][i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, 
						"neq"+to_string(r)+"_"+to_string(i)+"_"+to_string(j)+"_"+to_string(k));

			}
		}
	}

	//Sum of all the ones involved in the same round
	vector<GRBVar> sumneq_round(nbRound);
	for(uint r = 0; r < nbRound; r++)
		sumneq_round[r] = model.addVar(0.0,16*n,0.0,GRB_INTEGER, "sumneq_round"+to_string(r));
	//Sum of all neq vars
	GRBVar sumneq = model.addVar(0.0, nbRound*16*n, 0.0, GRB_INTEGER, "sumneq");

    model.update(); 
    GRBLinExpr expr_sumneq = 0;
    for(uint r = 0; r < nbRound; r++){
    	GRBLinExpr expr_sumneq_round = 0;
		for(uint i = 0; i < 4; i++){
			for(uint j = 0; j < 4; j++){
				for(uint k = 0; k < gamma-1; k++)
					expr_sumneq_round += neq[r][i][j][k];
				for(uint k = gamma; k < n-1; k++)
					expr_sumneq_round += neq[r][i][j][k];
			}
		}
		model.addConstr(sumneq_round[r] == expr_sumneq_round);
		expr_sumneq += sumneq_round[r];
	}
	model.addConstr(expr_sumneq == sumneq);
	model.update();

	//Constraints to link the neq vars
	addSalsaNeqConstr(model,nbRound,startingRound,gamma,x,z,neq);

	model.update();
	return model;
}

void
addModAddProbaSalsaQuarterRound(GRBModel & model,
								std::array<std::vector<GRBVar>,4> &  in,
								std::array<std::vector<GRBVar>,4> & out,
								std::vector<std::vector<GRBVar>> & z,
								std::vector<GRBVar> & wt,
								uint const gamma,
								std::string const & index){

	addModAddRXProbaConstr(model, in[0], in[3],z[0],wt[0],gamma,"_r0"+index);
	addModAddRXProbaConstr(model,out[1], in[0],z[1],wt[1],gamma,"_r1"+index);
	addModAddRXProbaConstr(model,out[2],out[1],z[2],wt[2],gamma,"_r2"+index);
	addModAddRXProbaConstr(model,out[3],out[2],z[3],wt[3],gamma,"_r3"+index);

}

std::vector<std::vector<std::vector<GRBVar>>>
addModAddProbaSalsa(GRBModel & model,
					uint const nbRound,
					uint const startingRound,
					uint const gamma,
					std::vector<std::vector<std::vector<GRBVar>>> & x,
					std::vector<std::vector<std::vector<std::vector<GRBVar>>>> & z){

	uint n = 32;

	//Weight variables
	vector<vector<vector<GRBVar>>> wt(nbRound, 
		   vector<vector<GRBVar>>    (4,
				  vector<GRBVar>	 (4)));
	for(uint r = 0; r < nbRound; r++){
		for(uint i = 0; i < 4; i++){
			for(uint j = 0; j < 4; j++){
				wt[r][i][j] = model.addVar(0,2*n,0,GRB_CONTINUOUS, 
					"wt"+to_string(r)+"_"+to_string(i)+"_"+to_string(j));
			}
		}
	}
	model.update();

	for(uint r = startingRound; r < startingRound+nbRound; r++){
		uint i = r-startingRound;
		auto &  xr = x[i];
		auto & xr1 = x[i+1];
		std::array<std::vector<GRBVar>,4>  in;
		std::array<std::vector<GRBVar>,4> out;

		if(r%2 == 0){
			//QR(x[ 0], x[ 4], x[ 8], x[12]);	// column 1
			 in[0] =  xr[0];  in[1] =  xr[4];  in[2] =  xr[8];  in[3] =  xr[12];
			out[0] = xr1[0]; out[1] = xr1[4]; out[2] = xr1[8]; out[3] = xr1[12];
			addModAddProbaSalsaQuarterRound(model,in,out,z[i][0],wt[i][0],gamma,"0"+to_string(i));

			//QR(x[ 5], x[ 9], x[13], x[ 1]);	// column 2
			 in[0] =  xr[5];  in[1] =  xr[9];  in[2] =  xr[13];  in[3] =  xr[1];
			out[0] = xr1[5]; out[1] = xr1[9]; out[2] = xr1[13]; out[3] = xr1[1];
			addModAddProbaSalsaQuarterRound(model,in,out,z[i][1],wt[i][1],gamma,"1"+to_string(i));

			//QR(x[10], x[14], x[ 2], x[ 6]);	// column 3
			 in[0] =  xr[10];  in[1] =  xr[14];  in[2] =  xr[2];  in[3] =  xr[6];
			out[0] = xr1[10]; out[1] = xr1[14]; out[2] = xr1[2]; out[3] = xr1[6];
			addModAddProbaSalsaQuarterRound(model,in,out,z[i][2],wt[i][2],gamma,"2"+to_string(i));

			//QR(x[15], x[ 3], x[ 7], x[11]);	// column 4
			 in[0] =  xr[15];  in[1] =  xr[3];  in[2] =  xr[7];  in[3] =  xr[11];
			out[0] = xr1[15]; out[1] = xr1[3]; out[2] = xr1[7]; out[3] = xr1[11];
			addModAddProbaSalsaQuarterRound(model,in,out,z[i][3],wt[i][3],gamma,"3"+to_string(i));
		}
		else{
			//QR(x[ 0], x[ 1], x[ 2], x[ 3]);	// row 1
			 in[0] =  xr[0];  in[1] =  xr[1];  in[2] =  xr[2];  in[3] =  xr[3];
			out[0] = xr1[0]; out[1] = xr1[1]; out[2] = xr1[2]; out[3] = xr1[3];
			addModAddProbaSalsaQuarterRound(model,in,out,z[i][0],wt[i][0],gamma,"0"+to_string(i));

			//QR(x[ 5], x[ 6], x[ 7], x[ 4]);	// row 2
			 in[0] =  xr[5];  in[1] =  xr[6];  in[2] =  xr[7];  in[3] =  xr[4];
			out[0] = xr1[5]; out[1] = xr1[6]; out[2] = xr1[7]; out[3] = xr1[4];
			addModAddProbaSalsaQuarterRound(model,in,out,z[i][1],wt[i][1],gamma,"1"+to_string(i));

			//QR(x[10], x[11], x[ 8], x[ 9]);	// row 3
			 in[0] =  xr[10];  in[1] =  xr[11];  in[2] =  xr[8];  in[3] =  xr[9];
			out[0] = xr1[10]; out[1] = xr1[11]; out[2] = xr1[8]; out[3] = xr1[9];
			addModAddProbaSalsaQuarterRound(model,in,out,z[i][2],wt[i][2],gamma,"2"+to_string(i));

			//QR(x[15], x[12], x[13], x[14]);	// row 4
			 in[0] =  xr[15];  in[1] =  xr[12];  in[2] =  xr[13];  in[3] =  xr[14];
			out[0] = xr1[15]; out[1] = xr1[12]; out[2] = xr1[13]; out[3] = xr1[14];
			addModAddProbaSalsaQuarterRound(model,in,out,z[i][3],wt[i][3],gamma,"3"+to_string(i));
		}
	}

	return wt;
}

void
addModAddVarsSalsa(std::array<std::vector<GRBVar>,4> &  in,
				   std::array<std::vector<GRBVar>,4> & out,
				   std::vector<std::vector<GRBVar>> & z,
				   std::vector<std::array<std::vector<GRBVar>,3>> & modAddVars,
				   uint & i){

	// addModAddRXProbaConstr(model, in[0], in[3],z[0],wt[0],gamma,"_r0"+index);
	modAddVars[i][0] = in[0]; modAddVars[i][1] = in[3]; modAddVars[i][2] = z[0]; i++;
	// addModAddRXProbaConstr(model,out[1], in[0],z[1],wt[1],gamma,"_r1"+index);
	modAddVars[i][0] =out[1]; modAddVars[i][1] = in[0]; modAddVars[i][2] = z[1]; i++;
	// addModAddRXProbaConstr(model,out[2],out[1],z[2],wt[2],gamma,"_r2"+index);
	modAddVars[i][0] =out[2]; modAddVars[i][1] =out[1]; modAddVars[i][2] = z[2]; i++;
	// addModAddRXProbaConstr(model,out[3],out[2],z[3],wt[3],gamma,"_r3"+index);
	modAddVars[i][0] =out[3]; modAddVars[i][1] =out[2]; modAddVars[i][2] = z[3]; i++;

}

std::vector<std::array<std::vector<GRBVar>,3>>
getModAddVarsSalsa(uint const nbRound,
				   uint const startingRound,
				   std::vector<std::vector<std::vector<GRBVar>>> & x,
				   std::vector<std::vector<std::vector<std::vector<GRBVar>>>> & z){
/*
Return a vector containing all triplets of variables that are involved in a modular addition, for the callback
*/

	vector<array<vector<GRBVar>,3>> modAddVars(16*nbRound);
	uint index = 0;
	for(uint r = startingRound; r < startingRound+nbRound; r++){
		uint i = r-startingRound;
		auto &  xr = x[i];
		auto & xr1 = x[i+1];
		std::array<std::vector<GRBVar>,4>  in;
		std::array<std::vector<GRBVar>,4> out;

		if(r%2 == 0){
			//QR(x[ 0], x[ 4], x[ 8], x[12]);	// column 1
			 in[0] =  xr[0];  in[1] =  xr[4];  in[2] =  xr[8];  in[3] =  xr[12];
			out[0] = xr1[0]; out[1] = xr1[4]; out[2] = xr1[8]; out[3] = xr1[12];
			addModAddVarsSalsa(in,out,z[i][0],modAddVars,index);

			//QR(x[ 5], x[ 9], x[13], x[ 1]);	// column 2
			 in[0] =  xr[5];  in[1] =  xr[9];  in[2] =  xr[13];  in[3] =  xr[1];
			out[0] = xr1[5]; out[1] = xr1[9]; out[2] = xr1[13]; out[3] = xr1[1];
			addModAddVarsSalsa(in,out,z[i][1],modAddVars,index);

			//QR(x[10], x[14], x[ 2], x[ 6]);	// column 3
			 in[0] =  xr[10];  in[1] =  xr[14];  in[2] =  xr[2];  in[3] =  xr[6];
			out[0] = xr1[10]; out[1] = xr1[14]; out[2] = xr1[2]; out[3] = xr1[6];
			addModAddVarsSalsa(in,out,z[i][2],modAddVars,index);

			//QR(x[15], x[ 3], x[ 7], x[11]);	// column 4
			 in[0] =  xr[15];  in[1] =  xr[3];  in[2] =  xr[7];  in[3] =  xr[11];
			out[0] = xr1[15]; out[1] = xr1[3]; out[2] = xr1[7]; out[3] = xr1[11];
			addModAddVarsSalsa(in,out,z[i][3],modAddVars,index);
		}
		else{
			//QR(x[ 0], x[ 1], x[ 2], x[ 3]);	// row 1
			 in[0] =  xr[0];  in[1] =  xr[1];  in[2] =  xr[2];  in[3] =  xr[3];
			out[0] = xr1[0]; out[1] = xr1[1]; out[2] = xr1[2]; out[3] = xr1[3];
			addModAddVarsSalsa(in,out,z[i][0],modAddVars,index);

			//QR(x[ 5], x[ 6], x[ 7], x[ 4]);	// row 2
			 in[0] =  xr[5];  in[1] =  xr[6];  in[2] =  xr[7];  in[3] =  xr[4];
			out[0] = xr1[5]; out[1] = xr1[6]; out[2] = xr1[7]; out[3] = xr1[4];
			addModAddVarsSalsa(in,out,z[i][1],modAddVars,index);

			//QR(x[10], x[11], x[ 8], x[ 9]);	// row 3
			 in[0] =  xr[10];  in[1] =  xr[11];  in[2] =  xr[8];  in[3] =  xr[9];
			out[0] = xr1[10]; out[1] = xr1[11]; out[2] = xr1[8]; out[3] = xr1[9];
			addModAddVarsSalsa(in,out,z[i][2],modAddVars,index);

			//QR(x[15], x[12], x[13], x[14]);	// row 4
			 in[0] =  xr[15];  in[1] =  xr[12];  in[2] =  xr[13];  in[3] =  xr[14];
			out[0] = xr1[15]; out[1] = xr1[12]; out[2] = xr1[13]; out[3] = xr1[14];
			addModAddVarsSalsa(in,out,z[i][3],modAddVars,index);
		}
	}
	return modAddVars;
}


std::tuple<double,uint, std::vector<std::array<uint64_t,16>>> 
searchBestTrailSalsa(uint const nbRound,
					 uint const startingRound,
					 uint const gamma,
					 GRBEnv & env,
					 bool const withExactProbability,
					 std::map<uint,uint> const & lowerBoundConsecutiveRounds,
					 std::vector<uint64_t> const & inputDiff,
					 std::vector<uint64_t> const & outputDiff){
	uint val_sumneq = 0;
	vector<array<uint64_t,16>> trail;
	auto knownTrail = make_pair(val_sumneq, trail);
	return searchBestTrailSalsa(nbRound,startingRound,gamma,env,withExactProbability,lowerBoundConsecutiveRounds,knownTrail,inputDiff,outputDiff);
}


std::tuple<double,uint, std::vector<std::array<uint64_t,16>>> 
searchBestTrailSalsa(uint const nbRound,
					 uint const startingRound,
					 uint const gamma,
					 GRBEnv & env,
					 bool const withExactProbability,
					 std::map<uint,uint> const & lowerBoundConsecutiveRounds,
					 std::pair<uint, std::vector<std::array<uint64_t,16>>> const & knownTrail,
					 std::vector<uint64_t> const & inputDiff,
					 std::vector<uint64_t> const & outputDiff){
/*
Search for the best RX-diff trail with RX rotation gamma
withExactProbability allows to control if we add probability constraints and optimization. False means we just optimize on the number of neq variables
lowerBoundConsecutiveRounds provides lower bounds on the number of Neq vars for consecutive rounds
e.g. if lowerBoundConsecutiveRounds[2] = 10, then we add constraints so that for any 2 consecutive rounds, there are at least 10 neq vars set to 1
knownTrail can be provided as a pair (W,trail) where W is a bound on the number of neq variables and trail is a starting solution for the model
inputDiff and outputDiff can be left empty to allow any input/output differential, or given to fix the input/output differential respectively
*/

	uint n = 32;
	auto model = getSalsaModel(nbRound,startingRound,gamma,env);

	vector<vector<vector<GRBVar>>> x(nbRound+1);
	vector<vector<vector<vector<GRBVar>>>> z(nbRound,vector<vector<vector<GRBVar>>>(4));
	for(uint r = 0; r < nbRound+1; r++)
		x[r] = get2DArrayVar(model, 16, n, "dx"+to_string(r)+"_");
	for(uint r = 0; r < nbRound; r++){
		for(uint i = 0; i < 4; i++)
			z[r][i] = get2DArrayVar(model, 4, n, "dz"+to_string(r)+"_"+to_string(i)+"_");
	}

	//--- tmp change ---
	// cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	// cout << "!!!!! Custom model, not the actual one, check \"Tmp change\" flag !!!" << endl;
	// cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	// //This is just for experimentation
	// GRBLinExpr expr_nonzero = 0;
	// for(uint i = 0; i < 16; i++){
	// 	for(uint j = 0; j < gamma-1; j++)
	// 		expr_nonzero += x[0][i][j];
	// 	for(uint j = gamma; j < n; j++)
	// 		expr_nonzero += x[0][i][j];
	// }
	// model.addConstr(expr_nonzero >= 1);

	// for(uint i = 0; i < 16; i++)
	// 	addValueToBinVectorConstr(model, x[nbRound][i], 0);
	//-------------------------------------------------------

	// --- Temporary solution for constants ---
	//Values of the constant over the diagonal
	array<uint32_t,4> salsaCst = {0x61707865, 0x3320646e, 0x79622d32, 0x6b206574};
	array<uint32_t,4> salsaCstRXDiff;
	for(uint i = 0; i < 4; i++)
		salsaCstRXDiff[i] = CSHL(salsaCst[i],gamma,n)^salsaCst[i];
	//Diagonal is x[0,5,10,15]
	array<uint,4> indexDiag = {0,5,10,15};
	for(uint i = 0; i < 4; i++)
		addValueToBinVectorConstr(model, x[0][indexDiag[i]], salsaCstRXDiff[i]);
	// ----------------------------------------


	//Input and output constraints (if any)
	if(inputDiff.size() > 0){
		for(uint i = 0; i < 16; i++)
			addValueToBinVectorConstr(model, x[0][i], inputDiff[i]);
	}
	if(outputDiff.size() > 0){
		for(uint i = 0; i < 16; i++)
			addValueToBinVectorConstr(model, x[nbRound][i], outputDiff[i]);
	}

	//Sum of all neq vars
    GRBVar sumneq = model.getVarByName("sumneq");

    //Prepare the callback
    auto modAddVars = getModAddVarsSalsa(nbRound,startingRound,x,z);
    vector<vector<GRBVar>> stateVars(nbRound+1, vector<GRBVar>(16*n));
    for(uint r = 0; r < nbRound+1; r++){
    	for(uint i = 0; i < 16; i++){
    		for(uint j = 0; j < n; j++){
    			stateVars[r][i*n+j] = x[r][i][j];
    		}
    	}
    }
    cout << "modAddVars.size() = " << modAddVars.size() << endl;

    GRBEnv tmpEnv = GRBEnv();
    tmpEnv.set(GRB_IntParam_OutputFlag, 0);
    SalsaCallback salsacb(n,gamma,modAddVars,stateVars,sumneq,nbRound,startingRound,tmpEnv);
    CustomCallback ccb(n,gamma,modAddVars,stateVars,sumneq);
	
    model.update();
    // model.setCallback(&salsacb);
    model.setCallback(&ccb);

    // if(withExactProbability) //Also check for a plaintext pair
	// 	model.setCallback(&salsacb);
	// else//Just intermediary print
	// 	model.setCallback(&ccb);

	//Check if we have a valid trail provided as input, if so, use it to bound the number of neq variables and provide a starting value for the trail variables
	auto const & trail = knownTrail.second;
	if(trail.size() > 0){
		// auto const bestBound = knownTrail.first;
		//Bound the number of notAllEq
		// model.addConstr(sumneq <= bestBound);

		//Provide MIP Start info from the the trail
		for(uint r = 0; r < nbRound+1; r++){
			auto const & trail_r = trail[r];
			for(uint i = 0; i < 16; i++){
				auto const & trail_ri = trail_r[i];
				for(uint j = 0; j < n; j++){
					uint val = (trail_ri >> j)&1;
					x[r][i][j].set(GRB_DoubleAttr_Start, val);
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
    	auto wt = addModAddProbaSalsa(model,nbRound,startingRound,gamma,x,z);
    	GRBLinExpr obj = 0;
		for(uint r = 0; r < nbRound; r++){
			for(uint i = 0; i < 4; i++){
				for(uint j = 0; j < 4; j++)
					obj += wt[r][i][j];
			}
		}
		//Objective: minimize probability
		model.setObjective(obj, GRB_MINIMIZE);
    }
    else{
    	//Objective: minimize the number of notAllEq
		GRBLinExpr obj = sumneq;
		model.setObjective(obj, GRB_MINIMIZE);
    }

    model.update();
    cout << "optimize ready" << endl; 
	model.optimize();
	if(model.get(GRB_IntAttr_SolCount) > 0){

		vector<vector<uint64_t>> valx(nbRound+1, vector<uint64_t>(16));
		vector<uint64_t> alpha(modAddVars.size());
		vector<uint64_t>  beta(modAddVars.size());
		vector<uint64_t> delta(modAddVars.size());
		vector<double> modAddLog(modAddVars.size());
		double sumLog = 0;

		for(uint r = 0; r < nbRound+1; r++){
			for(uint i = 0; i < 16; i++)
				valx[r][i] = getUintSolutionFromBinVector(x[r][i]);
		}
		for(uint i = 0; i < modAddVars.size(); i++){
			alpha[i] = getUintSolutionFromBinVector(modAddVars[i][0]);
			beta[i] = getUintSolutionFromBinVector(modAddVars[i][1]);
			delta[i] = getUintSolutionFromBinVector(modAddVars[i][2]);

			uint64_t countsol = getRXDiffCount(alpha[i],beta[i],delta[i],n,gamma);
			sumLog += (log2(countsol)-2*n);
			modAddLog[i] = log2(countsol)-2*n;
		}

		uint val_sumneq = uint(round(sumneq.get(GRB_DoubleAttr_X)));
		auto bestObj = model.get(GRB_DoubleAttr_ObjVal);
		auto bestLB = model.get(GRB_DoubleAttr_ObjBound);

		cout << "***********************************************************************" << endl;
		cout << "*** Salsa " << nbRound << " rounds with k = " << gamma << " ***" << endl;
		cout << "*** Best solution found with objective " << bestObj << " ***" << endl;
		if(model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
			cout << "Solution obtained after time-out" << endl;
		else
			cout << "Solution obtained after full optimization" << endl;

		cout << "Number of notAllEqual vars : " << val_sumneq << endl;
		cout << "Best lower bound on objective : " << bestLB << endl;
		cout << "Trail :" << endl;
		for(uint r = 0; r < nbRound+1; r++){
			cout << "x" << r << " : ";
			for(uint i = 0; i < 16; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << valx[r][i] << " ";
			cout << endl;
		}

		cout << "modAdd :" << endl;
		for(uint i = 0; i < alpha.size(); i++){
			cout << "0x" << setfill('0') << setw(n/4) << hex << alpha[i] << " + ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << beta[i] << " -> ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << delta[i] << dec;
			cout << " (2^" << modAddLog[i] << ")" << endl;
		}
		cout << "sumLog = " << sumLog << endl;
		cout << "modAddTrail = [";
		for(uint i = 0; i < alpha.size(); i++){
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
		cout << dec;
		cout << salsacb.nbSolFound << " possible trails examined in callback" << endl;
		cout << "***********************************************************************" << endl;
		cout << dec;

		vector<array<uint64_t,16>> trail(nbRound+1);
		for(uint r = 0; r < nbRound+1; r++){
			for(uint i = 0; i < 16; i++)
				trail[r][i] = valx[r][i];
		}

		return make_tuple(bestLB, val_sumneq, trail);
	}
	else{
		cout << "No solution" << endl;
		uint val_sumneq = 0;
		double bestLB = 0;
		vector<array<uint64_t,16>> trail;
		return make_tuple(bestLB, val_sumneq, trail);
	}
}

std::vector<GRBVar>
addSalsaOpValueConstr(GRBModel & model,
				 vector<GRBVar> & a,
				 vector<GRBVar> & b,
				 vector<GRBVar> & c,
				 vector<GRBVar> & out,
				 uint const rot,
				 std::string const & prefix){
/*
Add the constraints to modelize out = ((a + b) <<< rot) ^ c
Create temp variables named prefix+"_"+i for i in [0,31]
Return the vector of temp vars
*/

	//First, tmp = (a + b)
	auto tmp = genArrayBinVar(model,32,prefix);
	addModAddValueConstr(model,a,b,tmp,"carry"+prefix);
	//Rotate tmp
	vector<GRBVar> rotTmp(32);
	for(uint i = 0; i < 32; i++)
		rotTmp[i] = tmp[mod(i-rot,32)];
	//out = rotTmp ^ c
	for(uint i = 0; i < 32; i++)
		addXORConstr(model,rotTmp[i],c[i],out[i]);
	return tmp;
}

std::vector<std::vector<GRBVar>>
addSalsaQuarterRoundValueConstr(GRBModel & model,
						   std::array<std::vector<GRBVar>,4> & in,
						   std::array<std::vector<GRBVar>,4> & out,
						   std::string const & suffix){
/*
Add constraints between in and out to modelize the value propagation from in to out
in and out should be of size (4,32)
suffix is used for unique variable names
Return an array of 4 vectors of variables used to modelize the mod add
Quarter round is defined as 
input {a,b,c,d}
b ^= (a + d) <<< 7;
c ^= (b + a) <<< 9;
d ^= (c + b) <<< 13;
a ^= (d + c) <<< 18;
output {a,b,c,d}
*/
	std::vector<std::vector<GRBVar>> z(4);
	//b ^= (a + d) <<< 7 i.e. out[1] = ((in[0] + in[3]) <<< 7) ^ in[1]
	z[0] = addSalsaOpValueConstr(model,in[0],in[3],in[1],out[1],7,"z"+suffix+"_0_");

	//c ^= (b + a) <<< 9 i.e. out[2] = ((out[1] + in[0]) <<< 9) ^ in[2]
	z[1] = addSalsaOpValueConstr(model,out[1],in[0],in[2],out[2],9,"z"+suffix+"_1_");

	//d ^= (c + b) <<< 13 i.e. out[3] = ((out[2] + out[1]) <<< 13) ^ in[3]
	z[2] = addSalsaOpValueConstr(model,out[2],out[1],in[3],out[3],13,"z"+suffix+"_2_");

	//a ^= (d + c) <<< 18 i.e. out[0] = ((out[3] + out[2]) <<< 18) ^ in[0]
	z[3] = addSalsaOpValueConstr(model,out[3],out[2],in[0],out[0],18,"z"+suffix+"_3_");
	return z;
}

std::pair<std::vector<std::vector<std::vector<GRBVar>>>, 
		  std::vector<std::vector<std::vector<std::vector<GRBVar>>>>>
addSalsaValueConstr(GRBModel & model,
				  uint const nbRound,
				  uint const startingRound,
				  std::string const & suffix){
/*
Add constraints for RX-differential propagation for Salsa
Return the variables {x,z} used in the modelization
*/

	uint n = 32;
	vector<vector<vector<GRBVar>>> x(nbRound+1);
	vector<vector<vector<vector<GRBVar>>>> z(nbRound, vector<vector<vector<GRBVar>>>(4));
	for(uint r = 0; r < nbRound+1; r++)
		x[r] = genArray2DBinVar(model, 16, n, "x"+suffix+to_string(r)+"_");

	for(uint r = startingRound; r < startingRound+nbRound; r++){
		uint i = r-startingRound;
		auto &  xr = x[i];
		auto & xr1 = x[i+1];
		std::array<std::vector<GRBVar>,4>  in;
		std::array<std::vector<GRBVar>,4> out;

		if(r%2 == 0){
			//QR(x[ 0], x[ 4], x[ 8], x[12]);	// column 1
			 in[0] =  xr[0];  in[1] =  xr[4];  in[2] =  xr[8];  in[3] =  xr[12];
			out[0] = xr1[0]; out[1] = xr1[4]; out[2] = xr1[8]; out[3] = xr1[12];
			z[i][0] = addSalsaQuarterRoundValueConstr(model,in,out, suffix+to_string(i)+"0");

			//QR(x[ 5], x[ 9], x[13], x[ 1]);	// column 2
			 in[0] =  xr[5];  in[1] =  xr[9];  in[2] =  xr[13];  in[3] =  xr[1];
			out[0] = xr1[5]; out[1] = xr1[9]; out[2] = xr1[13]; out[3] = xr1[1];
			z[i][1] = addSalsaQuarterRoundValueConstr(model,in,out, suffix+to_string(i)+"1");

			//QR(x[10], x[14], x[ 2], x[ 6]);	// column 3
			 in[0] =  xr[10];  in[1] =  xr[14];  in[2] =  xr[2];  in[3] =  xr[6];
			out[0] = xr1[10]; out[1] = xr1[14]; out[2] = xr1[2]; out[3] = xr1[6];
			z[i][2] = addSalsaQuarterRoundValueConstr(model,in,out, suffix+to_string(i)+"2");

			//QR(x[15], x[ 3], x[ 7], x[11]);	// column 4
			 in[0] =  xr[15];  in[1] =  xr[3];  in[2] =  xr[7];  in[3] =  xr[11];
			out[0] = xr1[15]; out[1] = xr1[3]; out[2] = xr1[7]; out[3] = xr1[11];
			z[i][3] = addSalsaQuarterRoundValueConstr(model,in,out, suffix+to_string(i)+"3");
		}
		else{
			//QR(x[ 0], x[ 1], x[ 2], x[ 3]);	// row 1
			 in[0] =  xr[0];  in[1] =  xr[1];  in[2] =  xr[2];  in[3] =  xr[3];
			out[0] = xr1[0]; out[1] = xr1[1]; out[2] = xr1[2]; out[3] = xr1[3];
			z[i][0] = addSalsaQuarterRoundValueConstr(model,in,out, suffix+to_string(i)+"0");

			//QR(x[ 5], x[ 6], x[ 7], x[ 4]);	// row 2
			 in[0] =  xr[5];  in[1] =  xr[6];  in[2] =  xr[7];  in[3] =  xr[4];
			out[0] = xr1[5]; out[1] = xr1[6]; out[2] = xr1[7]; out[3] = xr1[4];
			z[i][1] = addSalsaQuarterRoundValueConstr(model,in,out, suffix+to_string(i)+"1");

			//QR(x[10], x[11], x[ 8], x[ 9]);	// row 3
			 in[0] =  xr[10];  in[1] =  xr[11];  in[2] =  xr[8];  in[3] =  xr[9];
			out[0] = xr1[10]; out[1] = xr1[11]; out[2] = xr1[8]; out[3] = xr1[9];
			z[i][2] = addSalsaQuarterRoundValueConstr(model,in,out, suffix+to_string(i)+"2");

			//QR(x[15], x[12], x[13], x[14]);	// row 4
			 in[0] =  xr[15];  in[1] =  xr[12];  in[2] =  xr[13];  in[3] =  xr[14];
			out[0] = xr1[15]; out[1] = xr1[12]; out[2] = xr1[13]; out[3] = xr1[14];
			z[i][3] = addSalsaQuarterRoundValueConstr(model,in,out, suffix+to_string(i)+"3");
		}
	}
	return make_pair(x,z);
}

GRBModel
getSalsaPlaintextPairModel(uint const nbRound, 
						   uint const startingRound,
						   uint const gamma,
						   GRBEnv & env,
						   std::vector<std::array<uint64_t,16>> const & trail){
//Create a MILP model for Salsa to find a plaintext pair matching a given trail

	uint n = 32;
	GRBModel model = GRBModel(env);

	//Constraints for the differential trail
	auto [x,z] = addSalsaRXDConstr(model,nbRound,startingRound,gamma);

	//Constraints for the plaintext values
	auto [x0,z0] = addSalsaValueConstr(model,nbRound,startingRound,"0");
	auto [x1,z1] = addSalsaValueConstr(model,nbRound,startingRound,"1");


	// --- Temporary solution for constants ---
	//Values of the constant over the diagonal
	array<uint32_t,4> salsaCst = {0x61707865, 0x3320646e, 0x79622d32, 0x6b206574};
	//Diagonal is x[0,5,10,15]
	array<uint,4> indexDiag = {0,5,10,15};
	for(uint i = 0; i < 4; i++){
		addValueToBinVectorConstr(model, x0[0][indexDiag[i]], salsaCst[i]);
		addValueToBinVectorConstr(model, x1[0][indexDiag[i]], salsaCst[i]);
	}
	// ----------------------------------------


	//Bind the value variables to the differential variables
	//Essentially, enforce dx = rot(x0) xor x1
	for(uint r = 0; r < nbRound+1; r++){
		for(uint i = 0; i < 16; i++)
			addRXDifferentialValueBinding(model,x0[r][i],x1[r][i],x[r][i],gamma);
	}
	for(uint r = 0; r < nbRound; r++){
		for(uint i = 0; i < 4; i++){
			for(uint j = 0; j < 4; j++)
				addRXDifferentialValueBinding(model,z0[r][i][j],z1[r][i][j],z[r][i][j],gamma);
		}
	}

	//Add the constraint values for the trail
	for(uint r = 0; r < nbRound+1; r++){
		for(uint i = 0; i < 16; i++)
			addValueToBinVectorConstr(model,x[r][i],trail[r][i]);
	}

	//Some parameters that should help solving
    //We're only interested in getting one solution
    model.set(GRB_IntParam_SolutionLimit, 1);
    //Focus on finding solutions rather than optimization
    model.set(GRB_IntParam_MIPFocus, 1);
    //Still give an arbitrary objective as it's supposed to help the solver a bit
    GRBLinExpr obj = 0;
    for(uint i = 0; i < 16; i++){
    	for(uint j = 0; j < n; j++)
    		obj += x0[0][i][j];
    }
    model.setObjective(obj, GRB_MINIMIZE);
    return model;
}

#define ROTL(a,b) (((a) << (b)) | ((a) >> (32 - (b))))

void QRSalsa(uint32_t & a,
			 uint32_t & b,
			 uint32_t & c,
			 uint32_t & d){
	b ^= ROTL(a + d, 7);
	c ^= ROTL(b + a, 9);
	d ^= ROTL(c + b,13);
	a ^= ROTL(d + c,18);
}

std::array<uint32_t,16>
salsaEncrypt(uint const nbRound,
			 uint const startingRound,
			 std::array<uint32_t,16> const & p){

	array<uint32_t,16> x(p);
	for(uint i = startingRound; i < startingRound+nbRound; i++){
		if(i%2 == 0){
			QRSalsa(x[ 0], x[ 4], x[ 8], x[12]);	// column 1
			QRSalsa(x[ 5], x[ 9], x[13], x[ 1]);	// column 2
			QRSalsa(x[10], x[14], x[ 2], x[ 6]);	// column 3
			QRSalsa(x[15], x[ 3], x[ 7], x[11]);	// column 4
		}
		else{
			QRSalsa(x[ 0], x[ 1], x[ 2], x[ 3]);	// row 1
			QRSalsa(x[ 5], x[ 6], x[ 7], x[ 4]);	// row 2
			QRSalsa(x[10], x[11], x[ 8], x[ 9]);	// row 3
			QRSalsa(x[15], x[12], x[13], x[14]);	// row 4
		}
	}
	return x;
}


bool
findPlaintextForTrailSalsa(uint const nbRound, 
						   uint const startingRound,
						   uint const gamma,
						   GRBEnv & env,
						   std::vector<std::array<uint64_t,16>> const & trail,
						   bool const printOutput){

	auto model = getSalsaPlaintextPairModel(nbRound,startingRound,gamma,env,trail);
	model.optimize();

	//Extract the solution, if any
	if(model.get(GRB_IntAttr_SolCount) == 0){
		if(printOutput){
			cout << "*************************************" << endl;
			cout << "*** No matching pair of plaintext ***" << endl;
			cout << "*************************************" << endl;
		}
		return false;
	}
	else{
		uint n = 32;
		auto x0 = get2DArrayVar(model,16,n,"x00_");
		auto x1 = get2DArrayVar(model,16,n,"x10_");

		array<uint32_t,16> valx0;
		array<uint32_t,16> valx1;
		for(uint i = 0; i < 16; i++){
			valx0[i] = getUintSolutionFromBinVector(x0[i]);
			valx1[i] = getUintSolutionFromBinVector(x1[i]);
		}

		//Compute input diff
		array<uint32_t,16> dx;
		for(uint i = 0; i < 16; i++)
			dx[i] = CSHL(valx0[i],gamma,n) ^ valx1[i];
		//Encrypt
		auto cx0 = salsaEncrypt(nbRound,startingRound,valx0);
		auto cx1 = salsaEncrypt(nbRound,startingRound,valx1);
		//Compute output diff
		array<uint32_t,16> cdx;
		for(uint i = 0; i < 16; i++)
			cdx[i] = CSHL(cx0[i],gamma,n) ^ cx1[i];

		if(printOutput){
			cout << "*******************************" << endl;
			cout << "*** Found pair of plaintext ***" << endl;
			cout << "x0 = ";
			for(uint i = 0; i < 16; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << valx0[i] << " ";
			cout << endl;
			cout << "x1 = ";
			for(uint i = 0; i < 16; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << valx1[i] << " ";
			cout << endl;

			cout << "Resulting input dx: ";
			for(uint i = 0; i < 16; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << dx[i] << " ";
			cout << endl;
			cout << "Expected input dx : ";
			for(uint i = 0; i < 16; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << trail[0][i] << " ";
			cout << endl;

			cout << "Resulting output dx: ";
			for(uint i = 0; i < 16; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << cdx[i] << " ";
			cout << endl;
			cout << "Expected output dx : ";
			for(uint i = 0; i < 16; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << trail[nbRound][i] << " ";
			cout << endl;
			cout << dec;


			for(uint i = 0; i < 16; i++){
				if(dx[i] != trail[0][i])
					cout << "Input differentials fail" << endl;
			}
			for(uint i = 0; i < 16; i++){
				if(cdx[i] != trail[nbRound][i])
					cout << "Output differentials fail" << endl;
			}
			cout << "*******************************" << endl;
		}

		return true;
	}
}






SalsaCallback::SalsaCallback(uint const wordSize,
							  uint const rotationOffset,
							  std::vector<std::array<std::vector<GRBVar>,3>> const & modAddVars,
							  std::vector<std::vector<GRBVar>> const & stateVars,
							  GRBVar & sumneq,
							  uint const nbRound,
							  uint const startingRound,
							  GRBEnv & env):
							  wordSize(wordSize),
							  rotationOffset(rotationOffset),
							  modAddVars(modAddVars),
							  stateVars(stateVars),
							  sumneq(sumneq),
							  nbRound(nbRound),
							  startingRound(startingRound),
							  env(env),
							  nbSolFound(0)
{}

void SalsaCallback::callback(){
try{
	if (where == GRB_CB_MIPSOL) { //If the solver found a solution

		nbSolFound++;
		if(nbSolFound%100 == 0){
			cout << nbSolFound << " callback checks so far" << endl;
		}
		//Get the value of the state for the trail
		std::vector<std::array<uint64_t,16>> val_state(nbRound+1);
		for(uint r = 0; r < nbRound+1; r++){
			for(uint i = 0; i < 16; i++){
				uint64_t wordVal = 0;
				for(uint64_t j = 0; j < wordSize; j++){
					uint tmp = uint(round(getSolution(stateVars[r][i*wordSize + j])));
					if(tmp)
						wordVal |= (1ULL << j);
				}
				val_state[r][i] = wordVal;
			}
		}

		//Before doing anything else, check if there is a matching plaintext pair following this trail
		if(findPlaintextForTrailSalsa(nbRound,startingRound,rotationOffset,env,val_state,false)){
			//Found a plaintext, so we assume it's ok

			uint nbModAdd = modAddVars.size();
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

			uint val_sumneq = uint(round(getSolution(sumneq)));

			//Print the trail
			cout << "=========================================================" << endl;
			cout << "======== Solution found with Objective : " << getDoubleInfo(GRB_CB_MIPSOL_OBJ) << " ========" << endl;
			cout << "Number of notAllEqual vars : " << val_sumneq << endl;
			cout << "Trail :" << endl;
			for(uint r = 0; r < nbRound+1; r++){
				cout << "Round " << r << " : ";
				for(auto const & word : val_state[r])
					cout << "0x" << setfill('0') << setw(wordSize/4) << hex << word << dec << " ";
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
		else{
			// cout << "Possible Trail but failed plaintext checking :" << endl;
			// for(uint r = 0; r < nbRound+1; r++){
			// 	cout << "Round " << r << " : ";
			// 	for(auto const & word : val_state[r])
			// 		cout << "0x" << setfill('0') << setw(wordSize/4) << hex << word << dec << " ";
			// 	cout << endl;
			// }
			//No plaintext found, so bad trail
			GRBLinExpr cutExpr = 0;
			for(uint r = 0; r < nbRound+1; r++){
				for(uint i = 0; i < 16; i++){
					for(uint64_t j = 0; j < wordSize; j++){
						uint tmp = uint(round(getSolution(stateVars[r][i*wordSize + j])));
						if(tmp)
							cutExpr += (1 - stateVars[r][i*wordSize + j]);
						else
							cutExpr += stateVars[r][i*wordSize + j];
					}
				}
			}
			addLazy(cutExpr >= 1);
		}
	}
} catch (GRBException e) {
cout << "Error number: " << e.getErrorCode() << endl;
cout << e.getMessage() << endl;
} catch (...) {
cout << "Error during callback" << endl;
}
}





std::vector<std::vector<GRBVar>>
addSalsaQuarterRoundConstr(GRBModel & model,
						   uint const gamma,
						   std::vector<std::vector<GRBVar>> & in,
						   std::vector<std::vector<GRBVar>> & out,
						   std::string const & suffix){
	array<vector<GRBVar>,4> ain;
	copy_n(in.begin(), 4, ain.begin());
	array<vector<GRBVar>,4> aout;
	copy_n(out.begin(), 4, aout.begin());
	return addSalsaQuarterRoundConstr(model,gamma,ain,aout,suffix);
}

void 
addSalsaQuarterRoundNeqConstr(GRBModel & model,
							  std::vector<std::vector<GRBVar>> & in,
							  std::vector<std::vector<GRBVar>> & out,
							  std::vector<std::vector<GRBVar>> & z,
							  std::vector<std::vector<GRBVar>> & neq,
							  uint const gamma){
	array<vector<GRBVar>,4> ain;
	copy_n(in.begin(), 4, ain.begin());
	array<vector<GRBVar>,4> aout;
	copy_n(out.begin(), 4, aout.begin());
	return addSalsaQuarterRoundNeqConstr(model,ain,aout,z,neq,gamma);
}


GRBModel
getSalsaQuarterRoundModel(uint const gamma,
						  GRBEnv & env){
//Create a MILP model for just one quarter round of Salsa

	uint n = 32;
	GRBModel model = GRBModel(env);
	vector<vector<vector<GRBVar>>> x(2);
	x[0] = genArray2DBinVar(model, 4, n, "dx0_");
	x[1] = genArray2DBinVar(model, 4, n, "dx1_");
	model.update();
	auto z = addSalsaQuarterRoundConstr(model,gamma,x[0],x[1],"");

	//Neq vars
	vector<vector<GRBVar>> neq(4, vector<GRBVar>(n-1));
	for(uint j = 0; j < 4; j++){
		for(uint k = 0; k < gamma-1; k++)
			neq[j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, 
				"neq"+to_string(j)+"_"+to_string(k));
		for(uint k = gamma; k < n-1; k++)
			neq[j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, 
				"neq"+to_string(j)+"_"+to_string(k));
	}

	//Sum of all neq vars
	GRBVar sumneq = model.addVar(0.0, 4*n, 0.0, GRB_INTEGER, "sumneq");
	GRBLinExpr expr_sumneq = 0;
	for(uint j = 0; j < 4; j++){
		for(uint k = 0; k < gamma-1; k++)
			expr_sumneq += neq[j][k];
		for(uint k = gamma; k < n-1; k++)
			expr_sumneq += neq[j][k];
	}
	model.addConstr(expr_sumneq == sumneq);
	model.update();

	//Constraints to link the neq vars
	addSalsaQuarterRoundNeqConstr(model,x[0],x[1],z,neq,gamma);

	model.update();
	return model;
}

void
addModAddProbaSalsaQuarterRound(GRBModel & model,
								std::vector<std::vector<GRBVar>> & in,
								std::vector<std::vector<GRBVar>> & out,
								std::vector<std::vector<GRBVar>> & z,
								std::vector<GRBVar> & wt,
								uint const gamma,
								std::string const & index){
	array<vector<GRBVar>,4> ain;
	copy_n(in.begin(), 4, ain.begin());
	array<vector<GRBVar>,4> aout;
	copy_n(out.begin(), 4, aout.begin());
	addModAddProbaSalsaQuarterRound(model,ain,aout,z,wt,gamma,index);
}

std::vector<GRBVar>
addModAddProbaSalsaQR(GRBModel & model,
					  uint const gamma,
					  std::vector<std::vector<std::vector<GRBVar>>> & x,
					  std::vector<std::vector<GRBVar>> & z){
	uint n = 32;
	//weight variables
	vector<GRBVar> wt(4);
	for(uint j = 0; j < 4; j++)
		wt[j] = model.addVar(0,2*n,0,GRB_CONTINUOUS,"wt"+to_string(j));
	model.update();
	addModAddProbaSalsaQuarterRound(model,x[0],x[1],z,wt,gamma,"");
	return wt;
}

void
addModAddVarsSalsa(std::vector<std::vector<GRBVar>> & in,
				   std::vector<std::vector<GRBVar>> & out,
				   std::vector<std::vector<GRBVar>> & z,
				   std::vector<std::array<std::vector<GRBVar>,3>> & modAddVars,
				   uint & i){
	array<vector<GRBVar>,4> ain;
	copy_n(in.begin(), 4, ain.begin());
	array<vector<GRBVar>,4> aout;
	copy_n(out.begin(), 4, aout.begin());
	addModAddVarsSalsa(ain,aout,z,modAddVars,i);
}


std::tuple<double,uint, std::vector<std::array<uint64_t,4>>> 
searchBestTrailSalsaQR(int const indexCst,
					   uint const gamma,
					   GRBEnv & env,
					   bool const withExactProbability){
	uint val_sumneq = 0;
	vector<array<uint64_t,4>> trail;
	auto knownTrail = make_pair(val_sumneq, trail);
	return searchBestTrailSalsaQR(indexCst,gamma,env,withExactProbability,knownTrail);
}

std::tuple<double,uint, std::vector<std::array<uint64_t,4>>> 
searchBestTrailSalsaQR(int const indexCst,
					   uint const gamma,
					   GRBEnv & env,
					   bool const withExactProbability,
					   std::pair<uint, std::vector<std::array<uint64_t,4>>> const & knownTrail){
//Search for the best RX-diff trail over one quarter round 
//indexCst = -1 means we don't use the RXD of constants as constraints
//otherwise, fix the word of index indexCst to its corresponding RXD

	uint n = 32;
	uint nbRound = 1;
	auto model = getSalsaQuarterRoundModel(gamma,env);

	//Variables
	vector<vector<vector<GRBVar>>> x(2);
	x[0] = get2DArrayVar(model, 4, n, "dx0_");
	x[1] = get2DArrayVar(model, 4, n, "dx1_");
	auto z = get2DArrayVar(model, 4, n, "dz_");
	GRBVar sumneq = model.getVarByName("sumneq");



	// //-- tmp change
	// cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	// cout << "!!!!! Custom model, not the actual one, check \"Tmp change\" flag !!!" << endl;
	// cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	// GRBLinExpr tmpexpr = 0;
	// for(uint i = 0; i < 32; i++)
	// 	tmpexpr += x[0][0][i];
	// model.addConstr(tmpexpr >= 10);


	// for(uint i = 0; i < 4; i++)
	// 	addValueToBinVectorConstr(model, x[1][i], 0);

	// auto cst = genArrayBinVar(model,n,"cst");
	// addRXDifferentialValueBinding(model,cst,cst,x[0][0],gamma);

	// GRBLinExpr tmpexprcst = 0;
	// for(uint i = 0; i < 32; i++)
	// 	tmpexprcst += cst[i];
	// model.addConstr(tmpexprcst >= 12);
	// model.addConstr(tmpexprcst <= 20);
	// // -- end tmp change

	//Constraints from the constants, if needed
	if(indexCst >= 0 && indexCst < 4){
		//Values of the constant over the diagonal
		array<uint32_t,4> salsaCst = {0x61707865, 0x3320646e, 0x79622d32, 0x6b206574};
		//Compute the RXD
		array<uint32_t,4> salsaCstRXDiff;
		for(uint i = 0; i < 4; i++)
			salsaCstRXDiff[i] = CSHL(salsaCst[i],gamma,n)^salsaCst[i];
		//Add the RXD constraints
		addValueToBinVectorConstr(model, x[0][indexCst], salsaCstRXDiff[indexCst]);
	}

	//Prepare the callback
	vector<array<vector<GRBVar>,3>> modAddVars(4);
	uint tmpIndex = 0;
	addModAddVarsSalsa(x[0],x[1],z,modAddVars,tmpIndex);
	vector<vector<GRBVar>> stateVars(2, vector<GRBVar>(4*n));
	for(uint r = 0; r < 2; r++){
    	for(uint i = 0; i < 4; i++){
    		for(uint j = 0; j < n; j++){
    			stateVars[r][i*n+j] = x[r][i][j];
    		}
    	}
    }
    GRBEnv tmpEnv = GRBEnv();
	tmpEnv.set(GRB_IntParam_OutputFlag, 0);
	SalsaQRCallback salsacb(n,gamma,modAddVars,stateVars,sumneq,nbRound,tmpEnv,indexCst);
	CustomCallback ccb(n,gamma,modAddVars,stateVars,sumneq);

	model.update();
	model.setCallback(&salsacb);
	// model.setCallback(&ccb);

	//Check if we have a valid trail provided as input, if so, use it to bound the number of neq variables and provide a starting value for the trail variables
	auto const & trail = knownTrail.second;
	if(trail.size() > 0){
		// auto const bestBound = knownTrail.first;
		//Bound the number of notAllEq
		// model.addConstr(sumneq <= bestBound);

		//Provide MIP Start info from the the trail
		for(uint r = 0; r < nbRound+1; r++){
			auto const & trail_r = trail[r];
			for(uint i = 0; i < 4; i++){
				auto const & trail_ri = trail_r[i];
				for(uint j = 0; j < n; j++){
					uint val = (trail_ri >> j)&1;
					x[r][i][j].set(GRB_DoubleAttr_Start, val);
				}
			}
		}
	}

	//Set the objective depending on withExactProbability
	if(withExactProbability){
		auto wt = addModAddProbaSalsaQR(model,gamma,x,z);
		GRBLinExpr obj = 0;
		for(uint j = 0; j < 4; j++)
			obj += wt[j];

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

		vector<vector<uint64_t>> valx(nbRound+1, vector<uint64_t>(4));
		vector<uint64_t> alpha(modAddVars.size());
		vector<uint64_t>  beta(modAddVars.size());
		vector<uint64_t> delta(modAddVars.size());
		vector<double> modAddLog(modAddVars.size());
		double sumLog = 0;

		for(uint r = 0; r < nbRound+1; r++){
			for(uint i = 0; i < 4; i++)
				valx[r][i] = getUintSolutionFromBinVector(x[r][i]);
		}
		for(uint i = 0; i < modAddVars.size(); i++){
			alpha[i] = getUintSolutionFromBinVector(modAddVars[i][0]);
			beta[i] = getUintSolutionFromBinVector(modAddVars[i][1]);
			delta[i] = getUintSolutionFromBinVector(modAddVars[i][2]);

			uint64_t countsol = getRXDiffCount(alpha[i],beta[i],delta[i],n,gamma);
			sumLog += (log2(countsol)-2*n);
			modAddLog[i] = log2(countsol)-2*n;
		}

		uint val_sumneq = uint(round(sumneq.get(GRB_DoubleAttr_X)));
		auto bestObj = model.get(GRB_DoubleAttr_ObjVal);
		auto bestLB = model.get(GRB_DoubleAttr_ObjBound);

		cout << "***********************************************************************" << endl;
		cout << "*** SalsaQR cst " << indexCst << " with k = " << gamma << " ***" << endl;
		cout << "*** Best solution found with objective " << bestObj << " ***" << endl;
		if(model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
			cout << "Solution obtained after time-out" << endl;
		else
			cout << "Solution obtained after full optimization" << endl;

		cout << "Number of notAllEqual vars : " << val_sumneq << endl;
		cout << "Best lower bound on objective : " << bestLB << endl;
		cout << "Trail :" << endl;
		for(uint r = 0; r < nbRound+1; r++){
			cout << "x" << r << " : ";
			for(uint i = 0; i < 4; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << valx[r][i] << " ";
			cout << endl;
		}

		cout << "modAdd :" << endl;
		for(uint i = 0; i < alpha.size(); i++){
			cout << "0x" << setfill('0') << setw(n/4) << hex << alpha[i] << " + ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << beta[i] << " -> ";
			cout << "0x" << setfill('0') << setw(n/4) << hex << delta[i] << dec;
			cout << " (2^" << modAddLog[i] << ")" << endl;
		}
		cout << "sumLog = " << sumLog << endl;
		cout << "modAddTrail = [";
		for(uint i = 0; i < alpha.size(); i++){
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
		cout << dec;
		cout << salsacb.nbSolFound << " possible trails examined in callback" << endl;

		if(!withExactProbability){
			cout << "flagHeuristic: cst " << indexCst << " k = " << gamma << " nbNotAllEq = " << val_sumneq << " bestLB = " << bestLB << " sumLog = " << sumLog << endl; 
		}
		else{
			cout << "flagProbability: cst " << indexCst << " k = " << gamma << " nbNotAllEq = " << val_sumneq << " bestLB = " << bestLB << " sumLog = " << sumLog << endl; 
		}


		cout << "***********************************************************************" << endl;

		vector<array<uint64_t,4>> trail(nbRound+1);
		for(uint r = 0; r < nbRound+1; r++){
			for(uint i = 0; i < 4; i++)
				trail[r][i] = valx[r][i];
		}

		return make_tuple(bestLB, val_sumneq, trail);
	}
	else{
		cout << "No solution" << endl;
		uint val_sumneq = 0;
		double bestLB = 0;
		vector<array<uint64_t,4>> trail;
		return make_tuple(bestLB, val_sumneq, trail);
	}
}

std::vector<std::vector<GRBVar>>
addSalsaQuarterRoundValueConstr(GRBModel & model,
						   std::vector<std::vector<GRBVar>> & in,
						   std::vector<std::vector<GRBVar>> & out,
						   std::string const & suffix){
	array<vector<GRBVar>,4> ain;
	copy_n(in.begin(), 4, ain.begin());
	array<vector<GRBVar>,4> aout;
	copy_n(out.begin(), 4, aout.begin());
	return addSalsaQuarterRoundValueConstr(model,ain,aout,suffix);
}


GRBModel
getSalsaQRPlaintextPairModel(int const indexCst,
							 uint const gamma,
							 GRBEnv & env,
							 std::vector<std::array<uint64_t,4>> const & trail){
	uint n = 32;
	uint nbRound = 1;

	//Constraints for the differential trail
	GRBModel model = GRBModel(env);
	vector<vector<vector<GRBVar>>> x(2);
	x[0] = genArray2DBinVar(model, 4, n, "dx0_");
	x[1] = genArray2DBinVar(model, 4, n, "dx1_");
	model.update();
	auto z = addSalsaQuarterRoundConstr(model,gamma,x[0],x[1],"");

	//Constraints for the plaintext values
	vector<vector<vector<GRBVar>>> x0(2);
	x0[0] = genArray2DBinVar(model, 4, n, "x00_");
	x0[1] = genArray2DBinVar(model, 4, n, "x01_");
	model.update();
	auto z0 = addSalsaQuarterRoundValueConstr(model,x0[0],x0[1],"0");
	vector<vector<vector<GRBVar>>> x1(2);
	x1[0] = genArray2DBinVar(model, 4, n, "x10_");
	x1[1] = genArray2DBinVar(model, 4, n, "x11_");
	model.update();
	auto z1 = addSalsaQuarterRoundValueConstr(model,x1[0],x1[1],"1");

	// //--- Tmp change
	// cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	// cout << "!!!!! Custom model, not the actual one, check \"Tmp change\" flag !!!" << endl;
	// cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	// for(uint i = 0; i < 32; i++)
	// 	model.addConstr(x0[0][0][i] == x1[0][0][i]);


	// GRBLinExpr tmpexprcst = 0;
	// for(uint i = 0; i < 32; i++)
	// 	tmpexprcst += x0[0][0][i];

	// model.addConstr(tmpexprcst >= 12);
	// model.addConstr(tmpexprcst <= 20);
	// //-------

	//Constraints from the constants, if needed
	if(indexCst >= 0 && indexCst < 4){
		//Values of the constant over the diagonal
		array<uint32_t,4> salsaCst = {0x61707865, 0x3320646e, 0x79622d32, 0x6b206574};
		//Add the value constraints
		addValueToBinVectorConstr(model, x0[0][indexCst], salsaCst[indexCst]);
		addValueToBinVectorConstr(model, x1[0][indexCst], salsaCst[indexCst]);
	}

	//Bind the value variables to the differential variables
	//Essentially, enforce dx = rot(x0) xor x1
	for(uint r = 0; r < nbRound+1; r++){
		for(uint i = 0; i < 4; i++)
			addRXDifferentialValueBinding(model,x0[r][i],x1[r][i],x[r][i],gamma);
	}
	for(uint j = 0; j < 4; j++)
		addRXDifferentialValueBinding(model,z0[j],z1[j],z[j],gamma);

	//Add the constraint values for the trail
	for(uint r = 0; r < nbRound+1; r++){
		for(uint i = 0; i < 4; i++)
			addValueToBinVectorConstr(model,x[r][i],trail[r][i]);
	}

	//Some parameters that should help solving
    //We're only interested in getting one solution
    model.set(GRB_IntParam_SolutionLimit, 1);
    //Focus on finding solutions rather than optimization
    model.set(GRB_IntParam_MIPFocus, 1);
    //Still give an arbitrary objective as it's supposed to help the solver a bit
    GRBLinExpr obj = 0;
    for(uint i = 0; i < 4; i++){
    	for(uint j = 0; j < n; j++)
    		obj += x0[0][i][j];
    }
    model.setObjective(obj, GRB_MINIMIZE);
    return model;
}

bool
findPlaintextForTrailSalsaQR(int const indexCst,
							 uint const gamma,
							 GRBEnv & env,
							 std::vector<std::array<uint64_t,4>> const & trail,
							 bool const printOutput){

	auto model = getSalsaQRPlaintextPairModel(indexCst,gamma,env,trail);
	model.optimize();
	// uint n = 32;
	uint nbRound = 1;

	//Extract the solution, if any
	if(model.get(GRB_IntAttr_SolCount) == 0){
		if(printOutput){
			cout << "*************************************" << endl;
			cout << "*** No matching pair of plaintext ***" << endl;
			cout << "*************************************" << endl;
		}
		return false;
	}
	else{
		uint n = 32;
		auto x0 = get2DArrayVar(model,4,n,"x00_");
		auto x1 = get2DArrayVar(model,4,n,"x10_");

		array<uint32_t,4> valx0;
		array<uint32_t,4> valx1;
		for(uint i = 0; i < 4; i++){
			valx0[i] = getUintSolutionFromBinVector(x0[i]);
			valx1[i] = getUintSolutionFromBinVector(x1[i]);
		}

		//Compute input diff
		array<uint32_t,4> dx;
		for(uint i = 0; i < 4; i++)
			dx[i] = CSHL(valx0[i],gamma,n) ^ valx1[i];
		//Encrypt
		array<uint32_t,4> cx0(valx0);
		array<uint32_t,4> cx1(valx1);
		QRSalsa(cx0[0],cx0[1],cx0[2],cx0[3]);
		QRSalsa(cx1[0],cx1[1],cx1[2],cx1[3]);
		//Compute output diff
		array<uint32_t,4> cdx;
		for(uint i = 0; i < 4; i++)
			cdx[i] = CSHL(cx0[i],gamma,n) ^ cx1[i];

		if(printOutput){
			cout << "*******************************" << endl;
			cout << "*** Found pair of plaintext ***" << endl;
			cout << "x0 = ";
			for(uint i = 0; i < 4; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << valx0[i] << " ";
			cout << endl;
			cout << "x1 = ";
			for(uint i = 0; i < 4; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << valx1[i] << " ";
			cout << endl;

			cout << "Resulting input dx: ";
			for(uint i = 0; i < 4; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << dx[i] << " ";
			cout << endl;
			cout << "Expected input dx : ";
			for(uint i = 0; i < 4; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << trail[0][i] << " ";
			cout << endl;

			cout << "Resulting output dx: ";
			for(uint i = 0; i < 4; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << cdx[i] << " ";
			cout << endl;
			cout << "Expected output dx : ";
			for(uint i = 0; i < 4; i++)
				cout << "0x" << setfill('0') << setw(n/4) << hex << trail[nbRound][i] << " ";
			cout << endl;
			cout << dec;


			for(uint i = 0; i < 4; i++){
				if(dx[i] != trail[0][i])
					cout << "Input differentials fail" << endl;
			}
			for(uint i = 0; i < 4; i++){
				if(cdx[i] != trail[nbRound][i])
					cout << "Output differentials fail" << endl;
			}
			cout << "*******************************" << endl;
		}

		return true;
	}
}








SalsaQRCallback::SalsaQRCallback(uint const wordSize,
							  uint const rotationOffset,
							  std::vector<std::array<std::vector<GRBVar>,3>> const & modAddVars,
							  std::vector<std::vector<GRBVar>> const & stateVars,
							  GRBVar & sumneq,
							  uint const nbRound,
							  GRBEnv & env,
							  int const indexCst):
							  wordSize(wordSize),
							  rotationOffset(rotationOffset),
							  modAddVars(modAddVars),
							  stateVars(stateVars),
							  sumneq(sumneq),
							  nbRound(nbRound),
							  env(env),
							  nbSolFound(0),
							  indexCst(indexCst)
{}

void SalsaQRCallback::callback(){
try{
	if (where == GRB_CB_MIPSOL) { //If the solver found a solution

		nbSolFound++;
		if(nbSolFound%100 == 0){
			cout << nbSolFound << " callback checks so far" << endl;
		}
		//Get the value of the state for the trail
		std::vector<std::array<uint64_t,4>> val_state(nbRound+1);
		for(uint r = 0; r < nbRound+1; r++){
			for(uint i = 0; i < 4; i++){
				uint64_t wordVal = 0;
				for(uint64_t j = 0; j < wordSize; j++){
					uint tmp = uint(round(getSolution(stateVars[r][i*wordSize + j])));
					if(tmp)
						wordVal |= (1ULL << j);
				}
				val_state[r][i] = wordVal;
			}
		}

		//Before doing anything else, check if there is a matching plaintext pair following this trail
		if(findPlaintextForTrailSalsaQR(indexCst,rotationOffset,env,val_state,false)){
			//Found a plaintext, so we assume it's ok

			uint nbModAdd = modAddVars.size();
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

			uint val_sumneq = uint(round(getSolution(sumneq)));

			//Print the trail
			cout << "=========================================================" << endl;
			cout << "======== Solution found with Objective : " << getDoubleInfo(GRB_CB_MIPSOL_OBJ) << " ========" << endl;
			cout << "Number of notAllEqual vars : " << val_sumneq << endl;
			cout << "Trail :" << endl;
			for(uint r = 0; r < nbRound+1; r++){
				cout << "Round " << r << " : ";
				for(auto const & word : val_state[r])
					cout << "0x" << setfill('0') << setw(wordSize/4) << hex << word << dec << " ";
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
		else{
			// cout << "Possible Trail but failed plaintext checking :" << endl;
			// for(uint r = 0; r < nbRound+1; r++){
			// 	cout << "Round " << r << " : ";
			// 	for(auto const & word : val_state[r])
			// 		cout << "0x" << setfill('0') << setw(wordSize/4) << hex << word << dec << " ";
			// 	cout << endl;
			// }
			//No plaintext found, so bad trail
			GRBLinExpr cutExpr = 0;
			for(uint r = 0; r < nbRound+1; r++){
				for(uint i = 0; i < 4; i++){
					for(uint64_t j = 0; j < wordSize; j++){
						uint tmp = uint(round(getSolution(stateVars[r][i*wordSize + j])));
						if(tmp)
							cutExpr += (1 - stateVars[r][i*wordSize + j]);
						else
							cutExpr += stateVars[r][i*wordSize + j];
					}
				}
			}
			addLazy(cutExpr >= 1);
		}
	}
} catch (GRBException e) {
cout << "Error number: " << e.getErrorCode() << endl;
cout << e.getMessage() << endl;
} catch (...) {
cout << "Error during callback" << endl;
}
}


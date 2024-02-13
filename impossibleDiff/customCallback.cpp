#include "customCallback.hpp"

using namespace std;

typedef unsigned int uint;

CustomCallback::CustomCallback(std::vector<GRBVar> const & xkeyvar,
							   std::vector<unsigned int> const & xcutIndex):
							   keyvar(xkeyvar),
							   cutIndex(xcutIndex),
							   ctrKey(0),
							   foundKey()
{}

void CustomCallback::callback(){
	try{
	if (where == GRB_CB_MIPSOL) { //If the solver found a solution
		uint const nm = keyvar.size();
		//Extract the key value in this solution
		ctrKey++;
		vector<uint> bink(nm);
		for(uint i = 0; i < nm; i++){
			bink[i] = uint(round(getSolution(keyvar[i])));
		}
		// cout << "Found a solution with key    : ";
		// for(uint i = 0; i < nm; i++)
		// 	cout << bink[i];
		// cout << " (" << ctrKey << " keys so far)";
		// cout << endl;
		foundKey.emplace(bink);

		//Add a lazy constraint to remove any solution with this key value according to cutIndex
		GRBLinExpr cutExpr(0);
		vector<string> pattern(nm,"*");
		for(auto const i : cutIndex){
			if(bink[i] == 1) cutExpr += (1 - keyvar[i]);
			else cutExpr += keyvar[i];
			pattern[i] = to_string(bink[i]);
		}
		// cout << "Removed all keys with pattern: ";
		// for(auto const & tmp : pattern)
		// 	cout << tmp;
		// cout << " (" << ctrKey << " keys so far)";
		// cout << endl;
		
		addLazy(cutExpr >= 1);
	}

	} catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
	} catch (...) {
	cout << "Error during callback" << endl;
	}

}
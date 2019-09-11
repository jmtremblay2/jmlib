#include "model.hpp"
#include "random.hpp"
#include "arma_utils.hpp"

#include <algorithm>

using std::cout;
using std::endl;
using std::fill;

void rentix::Model::sampleFromPop(){
	jmlib::vdouble w = arma::toVector(weights);
	for(int i = 0; i < 10; ++i){
		jmlib::vint s = jmlib::sampleNoReplace(w,100);
	
		fill(useObs.begin(), useObs.end(), false);
		
		for(unsigned int j = 0; j < s.size(); ++j)
			useObs.at(s.at(j)) = true;
		this->solve(0);
	}
	
	return;
}

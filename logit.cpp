#include "logit.hpp"
#include "types.hpp"
//#include "jmlib_random.hpp"

#include <string>
#include <armadillo>

using std::string;
using std::cout;

rentix::Logit::Logit(string datafile, string model) :
        Model(model),
        Dlogit(datafile,model) {
	nobs = Dlogit.getNobs();
	nalt = Dlogit.getY().max();
	// Y takes values 0,1,...(nalt-1)
	nalt++;
	
	// number of parameters
	npar.push_back(Dlogit.getNvar());
	nparTotal = npar.at(0);

	// names of parameters
	names.push_back(jmlib::vstring(nparTotal));
	for(int i = 0; i < nparTotal; ++i)
		names.at(0).at(i) = Dlogit.getVarName(i);
	
	double bound = 25;
	jmlib::vvdouble low, up;
	up.push_back(jmlib::vdouble(npar.at(0),bound));
	low.push_back(jmlib::vdouble(npar.at(0),-bound));
	this->initSolverContainers(low,up,model);
	}

arma::colvec rentix::Logit::probas(const jmlib::vvdouble& x) const{
	const jmlib::vdouble& beta = x.at(0);
	arma::mat U = Dlogit.getUtilities(beta);
	
	// center the utilities at 0
	arma::colvec meanU = mean(U,1);
	for(unsigned int i = 0; i < U.n_cols; ++i)
		U.col(i) -= meanU;

	// compute probas
	U = arma::exp(U);
	arma::colvec chosenU(Dlogit.getNobs());
	arma::icolvec Y = Dlogit.getY();
	for(int obs = 0; obs < Dlogit.getNobs(); ++obs)
		chosenU(obs) = U(obs,Y(obs));

	arma::colvec sumU = arma::sum(U,1);
	
	return chosenU / sumU;
	}

/*	
void rentix::Logit::resampleY(const jmlib::vvdouble& x){
	const jmlib::vdouble& beta = x.at(0);
	arma::mat U = Dlogit.getUtilities(beta);
	for(int i = 0; i < U.n_rows; ++i)
		for(int j = 0; j < U.n_cols; ++j)
			U(i,j) += jmlib::rgumbel();
	
	for(int i = 0; i < nobs; ++i){
		// find the biggest utility;
		double biggest = U(i,0);
		int max_index = 0;
		for(int j = 1; j < U.n_cols; ++j)
			if(U(i,j) > biggest){
				biggest = U(i,j);
				max_index = j;
				}
		Dlogit.resetY(i,max_index);
		}
	return;
	}
*/

#include "model.hpp"
#include "oprobit.hpp"
#include "types.hpp"
#include "random.hpp"
#include "datareg.hpp"

#include <string>
#include <iostream>
#include <sstream>

using namespace std;

rentix::Oprobit::Oprobit(std::string datafile, std::string model):
			Model(model)
			//, Doprobit(datafile,model)
			{
	D = new rentix::DataReg(datafile, model);



	Doprobit = static_cast<DataReg*>(D);
	nalt = (int) (Doprobit->getY().max() + 0.5); // correction for float arithmetic
	// Y takes values 0,1,...,(nalt-1)
	nalt++;
	nobs = Doprobit->getNobs();

	// number of parameters
	npar.push_back(Doprobit->getNvar());
	npar.push_back(nalt-2);
	nparTotal = npar.at(0) + npar.at(1);

	// record the names of parameters to be estimated
	names.push_back(jmlib::vstring());
	for(int i = 0; i < npar.at(0); ++i)
		names.at(0).push_back(Doprobit->getVarName(i));

	names.push_back(jmlib::vstring());
	for(int i = 0; i < npar.at(1); ++i){
		std::stringstream ss;
		ss << i;
		names.at(1).push_back(std::string("alpha_") + ss.str());
		}
	double bound = 10000000;
	jmlib::vvdouble low(2), up(2);
	low.at(0) = jmlib::vdouble(npar.at(0), -bound);
	low.at(1) = jmlib::vdouble(npar.at(1), 0);
	up.at(0) = jmlib::vdouble(npar.at(0), bound);
	up.at(1) = jmlib::vdouble(npar.at(1), bound);
	// initialize the solver and containers
	this->initSolverContainers(low,up, model);
	}

arma::colvec rentix::Oprobit::probas(const jmlib::vvdouble& x) const {
	const jmlib::vdouble& beta = x.at(0);
	
	// recreate the gammas from the alphas (differences between gammas)
	jmlib::vdouble gamma(x.at(1));
	for(int i = 1; i < (nalt-2); ++i){
		gamma.at(i) += gamma.at(i-1);
		}
		
	arma::colvec Z = Doprobit->getXBeta(beta);
	arma::colvec probas(nobs);
	for(int obs = 0; obs < nobs; ++obs){
		// make sure Y is turned into an int without being rounded down
		int Y = (int) (Doprobit->getY(obs) + 0.5);
		// compute P(y)
		if(Y == 0)
			probas(obs) = jmlib::pnorm(-Z(obs));
		else if(Y == 1)
			probas(obs) = jmlib::pnorm(gamma.at(0) -Z(obs)) - jmlib::pnorm(-Z(obs));
		else if(1 < Y && Y < (nalt-1)){
			probas(obs) = jmlib::pnorm(gamma.at(Y-1) - Z(obs)) - jmlib::pnorm(gamma.at(Y-2) - Z(obs));
			}
		else // y == nalt - 1 : the highest alternative
			probas(obs) = 1 - jmlib::pnorm(gamma.at(nalt-3) - Z(obs));
		}
	return probas;
	}

// solve starts with 0 for the betas and 1 for the alphas
void rentix::Oprobit::solve(double def){
	jmlib::vvdouble start;
	// beta starts at 0
	start.push_back(jmlib::vdouble(npar.at(0)));
	for(int i = 0; i < npar.at(0); ++i)
		start.at(0).at(i) = def;

	// gamma starts at 1,2,...
	start.push_back(jmlib::vdouble(npar.at(1)));
	for(int i = 0; i < npar.at(1); ++i)
		start.at(1).at(i) = 1;

	Model::solve(start);
	return;
	}

void rentix::Oprobit::crap(){
/*	std::cout << "hello\n";	
	for(int i = 0; i < npar.size(); ++i)
		std::cout << npar.at(i) << std::endl;
	std::cout << "nalt : " << nalt << std::endl;*/
}

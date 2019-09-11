#include "dc2.hpp"

#include <iostream>
#include <string>
#include <armadillo>
#include <cmath>

#include "types.hpp"
#include "random.hpp"
#include "arma_utils.hpp"
#include "stdutils.hpp"
#include "utils.hpp"

using std::string;
using namespace std;

rentix::DC2::DC2(const std::string& datafile, const std::string& spec): Model(spec)
			, Dreg(datafile,spec,"regression","yreg"), Dop(datafile,spec,"op","yop"){
	nalt = ((int)(Dop.getY().max() + 0.5)) + 1;
	nparDisc = Dop.getNvar();
	nparCont = Dreg.getNvar();
	nobs = Dop.getNobs();

	npar.push_back(nparCont);
	npar.push_back(1);
	npar.push_back(nparDisc);
	npar.push_back(nalt-2);
	npar.push_back( 1);
	nparTotal = nparCont + 1 + nparDisc + (nalt-2) + 1;

	// names of parameters
	names = jmlib::vvstring(5);
	names.at(0) = jmlib::vstring(nparCont);
	for(int i = 0; i < nparCont; ++i)
		names.at(0).at(i) = Dreg.getVarName(i);
		
	names.at(1) = jmlib::vstring(1);
	names.at(1).at(0) = "sigma";
	
	names.at(2) = jmlib::vstring(nparDisc);
	for(int i = 0; i < nparCont; ++i)
		names.at(2).at(i) = Dop.getVarName(i);
	
	names.at(3) = jmlib::vstring(nalt - 2);
	for(int i = 0; i < nalt-2; ++i)
		names.at(3).at(i) = "alpha";
	
	names.at(4) = jmlib::vstring(1);
	names.at(4).at(0) = "rho";	
	
	
	// starting value
	start = jmlib::vvdouble(5);
	start.at(0) = jmlib::vdouble(nparCont);
	start.at(1) = jmlib::vdouble(1);
	start.at(2) = jmlib::vdouble(nparDisc);
	start.at(3) = jmlib::vdouble(nalt - 2);
	start.at(4) = jmlib::vdouble(1);
	fill(start.at(0).begin(), start.at(0).end(), 0);
	start.at(1).at(0) = 10;
	fill(start.at(2).begin(), start.at(2).end(), 0);
	fill(start.at(3).begin(), start.at(3).end(), 1);
	start.at(4).at(0) = 0;
	
	// bounds
	jmlib::vvdouble low(5), up(5);
	low.at(0) = up.at(0) = jmlib::vdouble(nparCont);
	low.at(1) = up.at(1) = jmlib::vdouble(1);
	low.at(2) = up.at(2) = jmlib::vdouble(nparDisc);
	low.at(3) = up.at(3) = jmlib::vdouble(nalt - 2);
	low.at(4) = up.at(4) = jmlib::vdouble(1);

	fill(low.at(0).begin(), low.at(0).end(), -1000000);
	low.at(1).at(0) = 0;
	fill(low.at(2).begin(), low.at(2).end(), -1000000);
	fill(low.at(3).begin(), low.at(3).end(), 0);
	low.at(4).at(0) = -1;
		
	fill(up.at(0).begin(), up.at(0).end(), 1000000);
	up.at(1).at(0) = 1000000;
	fill(up.at(2).begin(), up.at(2).end(), 1000000);
	fill(up.at(3).begin(), up.at(3).end(), 1000000);
	up.at(4).at(0) = 1;
	
	this->initSolverContainers(low,up,spec);
}

arma::colvec rentix::DC2::probas(const jmlib::vvdouble& x) const{
	arma::vec res;
	arma::vec regll = LLreg(x, res);
	arma::vec opll = LLop(x, res);
	
	arma::vec p(nobs);
	for(int i = 0; i < nobs; ++i)
		p(i) = regll(i) * opll(i);
	//return regll * opll;
	//return regll;
	return p;
}

arma::vec rentix::DC2::LLreg(const jmlib::vvdouble& x, arma::vec& e) const{
	arma::vec yhat = Dreg.getXBeta(x.at(0));
	e = Dreg.getY() - yhat;
	double sigma2 =  x.at(1).at(0);
	
	arma::vec l(nobs);
	for(int i = 0; i < nobs; ++i){
		l(i) = jmlib::dnorm(e(i), 0, sigma2);
		//std::cout << e(i) << "   " << l(i) << "   " << x.at(1).at(0) << "\n";
		}
	return l;
}

arma::vec rentix::DC2::LLop(const jmlib::vvdouble& x, const arma::vec& e) const{
	double rho = x.at(4).at(0);
	double sreg = std::sqrt(x.at(1).at(0));
	double scond = 1 - rho*rho;
	arma::vec mu = rho / sreg * e;
	
	// gamma
	jmlib::vdouble gamma = jmlib::vdouble(nalt + 1);
	gamma.at(0) = -100;
	gamma.at(1) = 0;
	for(int i = 0; i < nalt - 2; ++i)
		gamma.at(i + 2) = gamma.at(i + 1) + x.at(3).at(i);
	gamma.at(nalt) = 100;
	
	// conditional probit
	arma::vec op(nobs);
	arma::vec Z = Dop.getXBeta(x.at(2));
	for(int i = 0; i < nobs; ++i){
		int yop = (int) (Dop.getY(i) + 0.5);// + 1;
		double g1 = gamma.at(yop);
		double g2 = gamma.at(yop + 1);
		op(i) = jmlib::pnorm(g2 - Z(i), mu(i), scond) -  jmlib::pnorm(g1 - Z(i), mu(i), scond);
		//cout << yop << " | " << op(i) << " | " << i << "\n";
		}
	
	return op;
}

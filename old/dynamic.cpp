#include "dynamic.hpp"
#include <jm/types.hpp>
#include <jm/jmlib.hpp>
#include <jm/arma_utils.hpp>
#include <jm/constants.hpp>

#include <armadillo>
#include <iostream>

using std::cout;
using std::endl;

rentix::Dynamic::Dynamic(std::string model): Model(model){
	std::string fname = jmlib::getTokens("data",model).at(0).at(0);
	
	std::stringstream ss;
	ss << fname << ".choices.txt";
	std::string choiceName = ss.str();

	nTime = jmlib::countColumns(choiceName);
	nObs = jmlib::countLines(choiceName);
	nAlt = jmlib::getTokens("specific",model).size();
	nobs = nObs;

	choices = arma::toImat(jmlib::readMat(choiceName));

	
	D = std::vector<DataUt>(nTime + 2);
	for(int t = 0; t < nTime + 2; ++t){
		std::stringstream ss;
		ss << fname << "." << t << ".txt";
		D.at(t) = DataUt(ss.str(), model);
		//cout << jmlib::readMat(ss.str());
	}
	keeps = arma::imat(nObs, nTime);
	for(int obs = 0; obs < nObs; ++obs)
		for(int time = 0; time < nTime; ++time)
			keeps(obs,time) = choices(obs,time) == 0 ? 1 : 0;
	// get of 0 coded choices, now keep will have a -1 code

	npar = jmlib::vint(1);
	npar.at(0) = D.at(0).getNvar();
	nparTotal = npar.at(0);

	names.push_back(jmlib::vstring(npar.at(0)));
	for(int i = 0; i < npar.at(0); ++i)
		names.at(0).at(i) = D.at(0).getVarName(i);

	double bound = 25;
	jmlib::vvdouble up, low;
	up.push_back(jmlib::vdouble(npar.at(0),bound));
	low.push_back(jmlib::vdouble(npar.at(0),-bound));
	initSolverContainers(low,up,model);
}

void rentix::Dynamic::cinzia() const{

	jmlib::vvdouble x(1);
	x.at(0) = jmlib::vdouble(5,2);
	x.at(0).at(0) = -2;


	cout << "\n\nUTILITIES\n\n";

	
	for(int i = 0; i < nTime + 2; ++i)
		cout << "utilities for i = " << i << "\n"
				<< D.at(i).getUtilities(x.at(0));
	
	cout << "\n\nv_it\n\n";

	cout << getV(x);

	cout << "\n\nReservation\n\n";

	cout << getReservation(x);

	return;

	cout << "\n\nMode\n\n";

	cout << getMode(x);

	cout << "\n\nprobability of keeping\n\n";

	cout << getPKeep(x);

	cout << "\n\nprobability of alternatives\n\n";

	cout << getPAlts(x);

	cout << "\n\nprobability of decision\n\n";

	cout << getPDecision(x);

	cout << "\n\nproba for invidivuals\n\n";

	cout << probas(x);

	cout << "\n\none period payoff\n\n";

	cout << getOnePeriodPayoff(x);
}

arma::colvec rentix::Dynamic::probas(const jmlib::vvdouble& x) const {
	arma::mat PDecision = getPDecision(x);
	return prod(PDecision, 1);
}

arma::mat rentix::Dynamic::getPDecision(const jmlib::vvdouble& x) const {
	arma::mat P(nObs, nTime);
	arma::mat PKeep = getPKeep(x);
	arma::cube PAlts = getPAlts(x);

	for(int i = 0; i < nObs; ++i)
		for(int t = 0; t < nTime; ++t)//outMarket.at(i); ++t)
			if(1 == keeps(i,t))
				P(i,t) = PKeep(i,t);
			else
				P(i,t) = (1 - PKeep(i,t)) * PAlts(i, choices(i,t), t);
	return P;
}

arma::mat rentix::Dynamic::getPKeep(const jmlib::vvdouble& x) const{
	arma::mat Reservation = getReservation(x);
	arma::mat Mode = getMode(x);
	return exp( - exp( - (Reservation - Mode)));
}

arma::cube rentix::Dynamic::getPAlts(const jmlib::vvdouble& x) const{
	arma::cube PAlts(nObs, nAlt, nTime);
	for(int t = 0; t < nTime; ++t){
		arma::mat ExpUt = exp(D.at(t).getUtilities(x.at(0)));
		arma::colvec sumExpU = sum(ExpUt, 1);
		for(int alt = 0; alt < nAlt; ++alt)
			PAlts.slice(t).col(alt) = ExpUt.col(alt) / sumExpU;
		}
	return PAlts;
}

arma::mat rentix::Dynamic::getReservation(const jmlib::vvdouble& x) const{
	arma::mat onePerPayoff = getOnePeriodPayoff(x);
	arma::mat V = getV(x);
	arma::mat W(nObs, nTime);
	//cout << W;
	for(int i = 0; i < nObs; ++i)
		for(int t = 0; t < nTime; ++t){
			//W(i,t) = std::max(V(i,t), onePerPayoff(i,t) + std::max(V(i,t+1), onePerPayoff(i, t+1)));
			W(i,t) = std::max(V(i,t+1), V(i, t+2));
			//W(i,t) = std::max(W(i,t), 0);
			if(W(i,t) < 0)
				W(i,t) = 0;
		}
	//cout << W;
	return W;
}

arma::mat rentix::Dynamic::getMode(const jmlib::vvdouble& x) const {
	//arma::mat Mode(nObs, nTime);
	//Mode.fill(0);

	/*
	// old empirical mode
	arma::mat V = getV(x);
	arma::mat Mean = mean(V,0);
	Mean = Mean.submat(0,0,0, nTime-1);

	arma::mat R(nObs, nTime);
	for(int obs = 0; obs < nObs; ++obs)
		R.row(obs) = Mean - jmlib::gammaEuler;
	return R;
	//std::cout << Mean;
	//std::cout << R.row(0);
	*/
	
	arma::mat R(nObs, nTime);
	for(int t = 0; t < nTime; ++t){
		//arma::mat expV = arma::exp(D.at(t).getUtilities(x.at(0)));
		//std::cout << V.row(0);
		//arma::colvec sumExpU = arma::sum(V,1);
		//R.col(t) = arma::log(sumExpU);
		R.col(t) = arma::log(arma::sum(arma::exp(D.at(t).getUtilities(x.at(0))),1));
	}

	//std::cout << R.row(0);
	//int q;
	//std::cin >> q;
	return R;
}

arma::mat rentix::Dynamic::getOnePeriodPayoff(const jmlib::vvdouble& x) const{
	arma::mat oneP(nObs, nTime + 1);
	oneP.fill(0);
	return oneP;
}

arma::mat rentix::Dynamic::getV(const jmlib::vvdouble& x) const{
	// we need to go until nTime (and not nTime -1) because
	// we need to look one time period ahead
	arma::mat V(nObs, nTime + 2);
	for(int t = 0; t < nTime + 2; ++t){
		V.col(t) = max(D.at(t).getUtilities(x.at(0)), 1);
	}
	return V;
}

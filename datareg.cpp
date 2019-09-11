#include "types.hpp"
#include "jmlib.hpp"

#include "datareg.hpp"
#include "data.hpp"

#include <iostream>
#include <string>
#include <armadillo>

using std::cout;
using std::endl;

rentix::DataReg::DataReg(std::string datafile, std::string model,
		std::string reg_token, std::string Y_token) : Data(datafile) {
	// will contain only one vstring that will contain only one string (the name of the 
	// dependant variable)
	jmlib::vvstring vY = jmlib::getTokens(Y_token,model);
	// will contain only one vstring that will contain all the regressors
	jmlib::vvstring vReg = jmlib::getTokens(reg_token,model);

	// record name of Y and its values
	nameY = vY.at(0).at(0);
	Y = getVar(nameY);

	// record variables associated with indexes
	// and map from string to indexes
	index = vReg.at(0);
	nvar = index.size();
	for(int i = 0; i < nvar; ++i)
		names[index.at(i)] = i;

	// fill the X matrix
	dataReg = arma::mat(nobs,nvar);
	for(int obs = 0; obs < nobs; ++obs)
		for(int var= 0; var < nvar; ++var)
			dataReg(obs,var) = (*this)(obs,index.at(var)); // (*this) is a Data and can use operator()
	}

arma::rowvec rentix::DataReg::getObs(int obs) const {
	return dataReg.row(obs);
	}

double rentix::DataReg::getVal(int obs, int j) const {
    return dataReg(obs,j);
    }

arma::colvec rentix::DataReg::getXBeta(const arma::colvec& beta) const{
	return dataReg * beta;
	}

arma::colvec rentix::DataReg::getXBeta(const std::vector<double>& beta) const{
	// must turn beta into an arma::colvec before computations
	arma::colvec beta_arma(beta.size());
	for(unsigned int i = 0; i < beta.size(); ++i)
		beta_arma(i) = beta.at(i);
	return dataReg * beta_arma;
	}

// formulas of the OLS
arma::colvec rentix::DataReg::getBetaOLS() const {
	return arma::inv(arma::trans(dataReg) * dataReg) * arma::trans(dataReg) * Y;
	}

arma::colvec rentix::DataReg::getBetaSD() const {
	arma::mat A = getSigma2() * arma::inv(arma::trans(dataReg) * dataReg);
	arma::mat S = A.diag();
	return arma::sqrt(S);
	}

double rentix::DataReg::getSigma2() const {
	arma::colvec residuals = (dataReg * getBetaOLS() - Y);
	double SS = arma::as_scalar(arma::sum(residuals % residuals));
	return SS / (nobs - 1);
	}

double rentix::DataReg::getSigma2SD() const {
	return std::sqrt(2 * getSigma2() * getSigma2() / (nobs - nvar));
	}

double rentix::DataReg::getY(int i) const {
	return Y(i);
	}

arma::colvec rentix::DataReg::getY() const {
	return Y;
	}

std::string rentix::DataReg::getVarName(int index) const {
	return this->index.at(index);
	}

int rentix::DataReg::getVarIndex(std::string varName) const {
	return names.find(varName)->second;
	}

int rentix::DataReg::getNvar() const {
	return nvar;
	}

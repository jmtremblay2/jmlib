#include "dataut.hpp"
#include "jmlib.hpp"
#include "types.hpp"

#include <armadillo>
#include <string>
#include <iostream>

rentix::DataUt::DataUt(std::string datafile, std::string model) 
			: Data(datafile){
	// variables that have a common coefficient across alternatives
	// and variables with specific coefficient for each alternative
	jmlib::vvstring vY = jmlib::getTokens("Ydisc",model);
	jmlib::vvstring vCommon = jmlib::getTokens("common",model);
	jmlib::vvstring vSpecific = jmlib::getTokens("specific",model);
	// record name of Y and its values
	nameY = vY.at(0).at(0);
	Y = arma::icolvec(nobs);
	for(int i = 0; i < nobs; ++i)
		// make sure that double -> int does not decrement Y
		Y(i) = ((int) (*this)(i,nameY) + 0.5); 

	// record name of variables in the order they appear
	// in the dataut matrix
	for(unsigned int i = 0; i < vCommon.size(); ++i)
		index.push_back(vCommon.at(i).at(0));
	
	for(unsigned int i = 0; i < vSpecific.size(); ++i)
		for(unsigned int j = 0; j < vSpecific.at(i).size(); ++j)
			index.push_back(vSpecific.at(i).at(j));

	// maps from var names to int (id)
	for(unsigned int i = 0; i < index.size(); ++i)
		names[index.at(i)] = i;

	// record name of dependant variable, number of parameters and alternatives
	nvar = index.size();
	nalt = vSpecific.size();

	// modify at your own risks
	dataut = arma::mat(nobs*nalt,nvar);
	dataut.fill(0);
	for(int obs = 0; obs < nobs; ++obs){
		int index = 0;
		// each coefficient that is common accross alternatives
		for(unsigned int var = 0; var < vCommon.size(); ++ var){
			for(int alt = 0; alt < nalt; ++alt)
				dataut(nalt*obs + alt, index) =
						(*this)(obs, vCommon.at(var).at(alt));
			// increment index only after all alternatives have been filled
			// for that position
			index++; 
			}
		// for each alternative we write into dataUt the variables
		// whose coefficient are alternative-specific
		for(int alt = 0; alt < nalt; ++alt)
			for(unsigned int var = 0; var < vSpecific.at(alt).size(); ++var){
				dataut(nalt*obs + alt, index) = (*this)(obs, vSpecific.at(alt).at(var));
				index++; // index must be incremented every time because coeffcients arent shared
				}
		} // observations
	} // constructor

arma::mat rentix::DataUt::getObs(int obs) const {
	return dataut.submat(nalt*obs, 0, nalt*(obs+1) -1, nvar-1);
	}

arma::mat rentix::DataUt::getUtilities(const jmlib::vdouble& beta) const {
	arma::colvec beta_vec(beta.size());
	for(unsigned int i = 0; i < beta.size(); ++i)
		beta_vec(i) = beta.at(i);

	//std::cout << "beta in get utilities : \n" << beta_vec;
	//std::cout << "x matrix : \n" << dataut << "***\n";
	arma::mat returnMat(nobs,nalt);
	arma::colvec U = dataut * beta_vec;
	// U is a big vector with all utilities for observation 1, all utilities for obs 2...
	// turn it into a matrix with line == obs, columns == alternatives
	for(int obs = 0; obs < nobs; ++obs)
		for(int alt = 0; alt < nalt; ++alt)
			returnMat(obs,alt) = U(obs*nalt + alt);

	return returnMat;
	}

int rentix::DataUt::getY(int i) const {
	return Y(i);
	}

arma::icolvec rentix::DataUt::getY() const {
	return Y;
	}

std::string rentix::DataUt::getVarName(int index) const {
	return this->index.at(index);
	}

int rentix::DataUt::getVarIndex(std::string varName) const {
	return names.find(varName)->second;
	}

int rentix::DataUt::getNvar() const {
	return nvar;
	}

int rentix::DataUt::getNalt() const {
	return nalt;
	}

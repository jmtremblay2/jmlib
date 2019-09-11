#ifndef DATA_REGRESSION
#define DATA_REGRESSION

#include "data.hpp"
#include "datareg.hpp"
#include "types.hpp"

#include <vector>
#include <string>
#include <armadillo>

namespace rentix{

class DataReg : public Data{
	public:
		// constructor needs a data file and a model specification
		DataReg(std::string datafile, std::string model, 
				std::string regName = "regression", std::string YName = "Ycont");
		DataReg () {};

		// get one observation (not a row from the full data set)
		arma::rowvec getObs(int obs) const;
		double getVal(int obs, int j) const;
		
		// get one Variable by index
		arma::colvec getVarRaw(int i) const;

		// compute X Beta from either a colum vector or a std::vector
		arma::colvec getXBeta(const arma::colvec& beta) const;
		arma::colvec getXBeta(const std::vector<double>& beta) const;
		
		// OLS stuff for a regression
		arma::colvec getBetaOLS() const;
		arma::colvec getBetaSD() const;
		double getSigma2() const;
		double getSigma2SD() const;

		// the dependant variable
		double getY(int i) const;
		arma::colvec getY() const;
		
		void resetY(int i, double newY) { Y(i) = newY; }
		
		// ID <==> name
		std::string getVarName(int index) const;
		int getVarIndex(std::string varName) const;

		// number of selected variables
		int getNvar() const;

	protected:
		// ID <==> name
		jmlib::mapStringInt names;
		jmlib::vstring index;
		// the X matrix
		arma::mat dataReg;
		// dependant variable
		arma::colvec Y;
		std::string nameY;
		
		int nvar;

	};

} // namespace rentix
#endif

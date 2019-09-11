#ifndef DATA_UTILITY
#define DATA_UTILITY

#include "data.hpp"
#include "types.hpp"

#include <string>
#include <armadillo>

namespace rentix{

class DataUt : public Data{
	public:
		// constructor needs a data file and a model specification
		DataUt(std::string datafile, std::string model);
		DataUt() {}
		
		// one observation is a matrix with (# alternative) lines
		// and (# coefficients) columns
		arma::mat getObs(int obs) const;
		
		// nobs x nalt matrix with line #i == utilities of obs #1
		arma::mat getUtilities(const jmlib::vdouble& beta) const;
		
		// Y is discrete, starts at 0
		int getY (int i) const;
		arma::icolvec getY() const;
		
		void resetY(int i, int newY){ Y(i) = newY; }
		// name <==> ID
		std::string getVarName(int index) const;
		int getVarIndex(std::string varName) const;
		
		int getNvar() const;
		int getNalt() const;
	protected:
		// name <==> ID
		jmlib::mapStringInt names;
		jmlib::vstring index;
		
		// stacked observations, dataut is (Nalt * Nobs  X  Ncoefficients)
		arma::mat dataut;
		// discrete Y starting at 0
		arma::icolvec Y;
		std::string nameY;
		int nvar, nalt;
		
	};

} // namespace rentix	

#endif

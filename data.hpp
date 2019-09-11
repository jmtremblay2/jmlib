#ifndef DATA_JM
#define DATA_JM

#include "types.hpp"

#include <string>
#include <armadillo>

namespace rentix{

/**
 * Data input class
 *
 * Data reads a data file separated by spaces.
 * the first line of the file must contain the 
 * names of the variables without quotes.
 * for example:
 * AGE  DISTANCE  SEX  EDUC
 * the remaining lines must contain a full observation
 * of these variables, coded in double, separated
 * by sapces, NAs are not available yet.
 */
class Data{
	public:
		/**
		 * Constructor takes a file passed as a string
		 *
		 * @param dataFileName the name of the file with path
		 */
		Data(std::string dataFileName);
		Data(){}

		/**
		 * print a few lines of the data set
		 *
		 * @param from the number of the first line we want
		 * @param to the number of the last line we want
		 */
		void print(int from, int to) const;

		// get one (full) observation or one varialble
		/**
		 * get one observation from the data set
		 *
		 * @param i the index of the observation we want
		 *
		 * @return a row vector with all the variables of the i-th obs
		 */
		arma::rowvec getRow(int i) const;

		/**
		 * get the i-th variable in the data set
		 *
		 * @param i the index of the variable that we want
		 *
		 * @return a column vector of all the obs for this variable
		 */
		arma::colvec getVarRaw(int i) const;

		/**
		 * get a varialbe specified by name in the data set
		 *
		 * @param varName the name of the variable as specified in data file
		 *
		 * @return a colum vector containing all the obs for this variable
		 */
		arma::colvec getVar(std::string varName) const;

		/**
		 * get a variable specified by index and obs id
		 *
		 * @param obsNum the id of the observation we want
		 * @param varNum the id of the variable
		 *
		 * @return the value as a double
		 */
		double getValRaw(int obsNum, int varNum) const;
		
		/**
		 * get a variable by name and obs is
		 *
		 * @param obsNum the id of the observation
		 * @param varName the name of the variable 
		 *
		 * @return the value as a double
		 */
		double operator()(int rowNum, const std::string& varName) const;
		//double& operator()(int rowNum, const std::string& varName);

		/**
		 * get the name of a variable from its index
		 *
		 * @param index the indes of the variable in the data set
		 * 
		 * @return the name of the index-th variable
		 */
		std::string getVarNameRaw(int index) const;

		/**
		 * get the index of a variable from its name
		 *
		 * @param varName the name of the variable
		 *
		 * @return the index of the variable
		 */
		int getVarIndexRaw(std::string varName) const;

		/**
		 * get the number of observations in the data set
		 *
		 * @return the number of observations
		 */
		int getNobs() const;

		/**
		 * get the number of variables in the data set
		 *
		 * @return the number of variables (columns) in the data set
		 */
		int getNvarRaw() const;

	protected:
		// correspondance between names and ID
		jmlib::mapStringInt namesRaw;
		jmlib::vstring indexRaw;
		
		// matrix to put the data in
		arma::mat data;
		
		// sizes
		int nobs;
		int nvarRaw;

	};

} // namespace rentix
#endif

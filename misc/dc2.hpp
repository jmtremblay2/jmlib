#ifndef DC2_HPP
#define DC2_HPP

#include <string>
#include <armadillo>
#include "types.hpp"
#include "model.hpp"
#include "datareg.hpp"

using std::string;

namespace rentix{

struct DC2param {
	arma::vec breg;
	double s;
	arma::vec bop;
	arma::vec gamma;
	double rho;
};

class DC2: public Model{
	public:
		DC2(){};
		DC2(const std::string& datafile, const std::string& spec);
		arma::colvec probas(const jmlib::vvdouble& x) const;
		~DC2() {}
		void crap(){
			probas(start);
		}
	protected:
		DC2param getParam(const arma::vec& b) const;
		
		arma::vec LLreg(const jmlib::vvdouble& x, arma::vec& e) const;
		arma::vec LLop(const jmlib::vvdouble& x, const arma::vec& e) const;
	
		int nparDisc, nparCont;
		int nalt;
		//arma::mat Xdisc, Xcont;
		arma::vec Ycont;
		arma::ivec Ydisc;
		DataReg Dreg, Dop;
};
}//namespace

#endif

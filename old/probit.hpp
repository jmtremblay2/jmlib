#ifndef PROBIT
#define PROBIT

#include "jmlib_types.hpp"
#include "model.hpp"
#include "dataut.hpp"

#include <armadillo>
#include <string>
#include <vector>

namespace rentix{

class Probit : public Model{
	public:
		Probit(std::string datafile, std::string model);
		double LL(const jmlib::vvdouble& x, jmlib::vvdouble& grad) const; // the one the solver uses
		double LL(const jmlib::vvdouble& x) const;

	protected:
        void computeSE(const jmlib::vdouble& x);
		int nsim;
		jmlib::vmat errors_ind;
		//arma::mat chol_var;
		DataUt Dprobit;
		int nalt;
	};

} // namespace rentix
#endif

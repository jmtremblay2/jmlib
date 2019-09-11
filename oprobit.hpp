#ifndef OPROBIT
#define OPROBIT

#include "types.hpp"
#include "model.hpp"

#include "datareg.hpp"

#include <string>

namespace rentix{

class Oprobit : public Model{
	public:
		Oprobit() {}
		Oprobit(std::string datafile, std::string model);
		arma::colvec probas(const jmlib::vvdouble& x) const;
		
		// override to start with an ordered gamma parameter
		void solve(double def = 0);

	protected:
		DataReg* Doprobit;
		//DataReg* Doprobit2;
		int nalt;
		void crap();	
	};

} // namespace rentix

#endif

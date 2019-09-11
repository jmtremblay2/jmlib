#ifndef LOGIT
#define LOGIT

#include "types.hpp"
#include "model.hpp"
#include "dataut.hpp"

#include <string>

namespace rentix{

class Logit : public Model{
	public:
		Logit() {}
		Logit(std::string datafile, std::string model);
		arma::colvec probas(const jmlib::vvdouble& x) const;

	protected:
		// void resampleY(const jmlib::vvdouble& x);
		DataUt Dlogit;
		int nalt;
	};

} // namespace rentix

#endif

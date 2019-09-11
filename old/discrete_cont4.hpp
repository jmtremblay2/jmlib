#ifndef DC4
#define DC4

#include "dataut.hpp"
#include "datareg.hpp"
#include "model.hpp"

#include <string>
#include <vector>
#include <armadillo>

namespace rentix{

class DC4 : public Model {

	public:

		DC4(std::string datafile, std::string model);
		arma::colvec probas(const jmlib::vvdouble& x) const;
		void solve(double def = 0);
		
	protected:
		// for constructor
		void fillParams(std::string model);
		void fillNpar();
		void fillNames();
		void fillBounds();
		void genErrors(int err_size);
		void fillStart();

		//void resampleY(const jmlib::vvdouble& x);
		//jmlib::vvdouble paramH2Comp(const jmlib::vvdouble& x) const;
		//jmlib::vvdouble paramComp2H(const jmlib::vvdouble& x) const;
		
		
		DataUt Ddisc;
		DataReg Dcont;

		jmlib::vvdouble start;
		//double defaultStartSD;
		
		int nalt;
		int nsim;
		double delta;
		jmlib::vmat errorsInd;
	};

} // namespace rentix

#endif

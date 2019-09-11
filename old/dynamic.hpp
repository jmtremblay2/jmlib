#ifndef DYNAMIC_HEADER
#define DYNAMIC_HEADER

#include <jm/model.hpp>
#include <jm/dataut.hpp>
#include <armadillo>

namespace rentix{

enum dynstate{
	IN_MARKET
	,OUT_MARKET
};

class Dynamic : public Model{
	public:
		Dynamic(std::string model);
		arma::colvec probas(const jmlib::vvdouble& x) const;

		void cinzia() const;

	protected:

		std::vector<rentix::DataUt> D;
		jmlib::vint outMarket;

		int nTime;
		int nAlt;
		int nObs;

		arma::imat choices;
		arma::imat keeps;

		arma::mat getPDecision(const jmlib::vvdouble& x) const;
		arma::mat getPKeep(const jmlib::vvdouble& x) const;
		arma::cube getPAlts(const jmlib::vvdouble& x) const;
		arma::mat getReservation(const jmlib::vvdouble& x) const;
		arma::mat getMode(const jmlib::vvdouble& x) const;
		arma::mat getOnePeriodPayoff(const jmlib::vvdouble& x) const;
		arma::mat getV(const jmlib::vvdouble& x) const;
};
} // namespace
#endif

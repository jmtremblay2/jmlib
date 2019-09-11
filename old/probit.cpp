#include "probit.hpp"
#include "jmlib.hpp"
#include "jmlib_types.hpp"
#include "jmlib_random.hpp"
#include "jmlib_gradient.hpp"

#include <string>
#include <iostream>
#include <sstream>
#include <string>

rentix::Probit::Probit(std::string datafile, std::string model) :
        Model(model),
        Dprobit(datafile,model){

	nalt = Dprobit.getY().max();
	// Y takes values 0,1,...(nalt-1)
	nalt++;

	// number of parameters
	npar.push_back(Dprobit.getNvar());
	// Ok now each line of the cholesky factor of the covariance matrix
	// execpt the first one whose element is set to 1
	for(int i = 2; i <= nalt; ++i)
		npar.push_back(i);
	nparTotal = 0;
	for(int i = 0; i < npar.size(); ++i)
		nparTotal += npar.at(i);

    // names of parameters
    names = jmlib::vvstring(npar.size());
    for(int i = 0; i < npar.size(); ++i)
        names.at(i) = jmlib::vstring(npar.at(i));

    for(int i = 0; i < npar.size(); ++i)
        names.at(0).at(i) = Dprobit.getVarName(i);

    for(int i = 2; i <= nalt; ++i)
        for(int j = 1; j <= i; ++j){
            std::stringstream ss;
            ss << "chol_{" << i << "," << j << "}";
            names.at(i-1).at(j-1) = ss.str();
            }
    // we start with bounds (that are totally arbitraty)
    double lb = -5, ub = 5;
    jmlib::vvdouble low(npar.size()), up(npar.size());
    for(int i = 0; i < npar.size(); ++i){
        low.at(i) = jmlib::vdouble(npar.at(i),lb);
        up.at(i) = jmlib::vdouble(npar.at(i),ub);
        }

	nsim = jmlib::getVal<int>("nsim", model, 1000);

	this->initSolverContainers(low,up);

    // initialize random terms
	for(int obs = 0; obs < Dprobit.getNobs(); ++obs){
		errors_ind.push_back(arma::mat(nalt,nsim));
		for(int i = 0; i < nalt; ++i)
			for(int j = 0; j < nsim; ++j)
				errors_ind.at(obs)(i,j) = jmlib::rnorm();
		}
    }

double rentix::Probit::LL(const jmlib::vvdouble& x) const {
    ncount_nograd++;
	arma::mat chol_var(Dprobit.getNalt(),Dprobit.getNalt());
	const jmlib::vdouble& beta(x.at(0));
    chol_var(0,0) = 1; // this is the parameter we set to identify the scale of our model
	for(int i = 1; i < Dprobit.getNalt(); ++i)
		for(int j = 0; j <= i; ++j)
			chol_var(i,j) = (x.at(i).at(j));

	arma::mat U = Dprobit.getUtilities(beta);
	double LL = 0;
	for(int obs = 0; obs < Dprobit.getNobs(); ++obs){
		arma::mat errors_dep = chol_var * errors_ind.at(obs);
		double n_success = 0;
		for(int sim = 0; sim < nsim; ++sim){
			int Y = Dprobit.getY(obs);
			bool add = true;
			for(int alt = 0; alt < Dprobit.getNalt(); ++alt)
				if( alt != Y &&  (U(obs,alt) + errors_dep(alt,sim))
						> (U(obs,Y) + errors_dep(Y,sim)))
					add = false;
			if(add)
				n_success++;
			}
		if(n_success > 0.5)
			LL += std::log( n_success / ((double) nsim) );
		else
			LL += std::log(0.000000000000001);
		}
	return LL;
	}

double rentix::Probit::LL(const jmlib::vvdouble& x, jmlib::vvdouble& grad) const{
    ncount_grad++;
	if(!grad.empty())
        grad = jmlib::gradient(this,x,0.01);
	return LL(x);
	}

void rentix::Probit::computeSE(const jmlib::vdouble& x){
    for(int i = 0; i < npar.size(); ++i)
        standardError.at(i) = jmlib::vdouble(npar.at(i), -1);
    }

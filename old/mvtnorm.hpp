#ifndef MVTNORM_HEADER
#define MVTNORM_HEADER

#include "f2c.h"

#include <armadillo>

namespace jmlib{

enum matorder_t {
	BYROW,
	BYCOLUMN
};

double pmvnorm(const arma::colvec& x, const arma::colvec& mu
		, const arma::mat& sigma);

#ifdef __cplusplus
extern "C" {
#endif

int mvtdst_(integer *n, integer *nu, doublereal *lower, 
	doublereal *upper, integer *infin, doublereal *correl, doublereal *
	delta, integer *maxpts, doublereal *abseps, doublereal *releps, 
	doublereal *error, doublereal *value, integer *inform__);
	
#ifdef __cplusplus
	}
#endif

}
#endif

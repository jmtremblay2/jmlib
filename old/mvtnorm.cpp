#include "mvtnorm.hpp"

#include <armadillo>
#include <cmath>

#include <jm/random.hpp>


/*
 * Multivariate Normal CDF
 *
 * @param x the point for which we want P(X < x)
 * @param mu the mean of the distribution
 * @param sigma the covariance of the distribution
 *
 * @return P(Normal(mu,Sigma) < x)
 */
double jmlib::pmvnorm(const arma::colvec& x, const arma::colvec& mu
		, const arma::mat& sigma){
	int n = sigma.n_rows;
	// X normalized to zero mean and variance 1
	arma::colvec xnorm(x);
	for(int i = 0; i < n; ++i)
		xnorm(i) = (x(i) - mu(i)) / std::sqrt(sigma(i,i));

	// calculate the elements above the diagonal of the 
	// correlation matrix that corresponds to sigma
	arma::mat corr = jmlib::cov2cor(sigma);
	arma::colvec corrVec(n*(n-1)/2);
	int index = 0;
	for(int i = 0; i < (n-1); ++i)
		for(int j = i + 1; j < n; ++ j)
			corrVec(index++) = corr(i,j);
	
	// code required to specify the external fortran routine
	// that we want P(Normal < x) (2 for each element)
	arma::icolvec infin(n);
	infin.fill(2);

	// lower bound is -infinity. -100 is close enough to that
	arma::colvec lower(n);
	lower.fill(-100);

	// mean of normalized distribution is 0, 
	arma::colvec zeroVec(n);
	zeroVec.fill(0);

	// parameters for the external routine. I am not sure what they 
	// all do. check out Alan Genz's website for more info. he is smarter
	// than me
	static int maxPoints = 2500;
	static int df = 0;
	static int inform = 0;
	static double abseps = 0.001;
	static double releps = 0;

	// estimated normal CDF along with estimated error.
	double error;
	double value;

	// external routine (in C but translated from fortran with f2c)
	mvtdst_(&n, &df, lower.memptr(), xnorm.memptr(), infin.memptr()
			, corrVec.memptr(), zeroVec.memptr(), &maxPoints, &abseps
			,&releps, &error, &value, &inform);

	return value;		
}

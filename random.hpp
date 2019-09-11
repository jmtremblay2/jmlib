#ifndef JMLIB_RANDOM
#define JMLIB_RANDOM

#include <boost/random/mersenne_twister.hpp>
#include <armadillo>
#include "types.hpp"

namespace jmlib{

/**
 * a random generating stream
 */
static boost::mt11213b mt; // stream

/**
 * random normal variable generation
 *
 * @param mu the mean of the normal
 * @param sigma2 the variable of the normal
 *
 * @return a random draw of the random variable
 */
double rnorm(double mu = 0, double sigma2 = 1);

/**
 * density of the random normal variable
 *
 * @param x the point to evaluate the density
 * @param mu the mean of the normal
 * @param sigma2 the variance of the normal
 *
 * @return the density of a N(mu,sigma2) evaluated at x
 */
double dnorm(double x, double mu = 0, double sigma2 = 1);

/**
 * CDF of the normal distribution
 *
 * @param x the point to evaluate
 * @param mu the mean of the normal
 * @param sigma2 the variable of the normal
 *
 * @return proba( N(mu,sigma2) < x)
 */
double pnorm(double x, double mu = 0, double sigma2 = 1);

/**
 * random draw of the gumbel distribution
 *
 * @return a draw of the random gumbel
 */
double rgumbel();

double runif();

vint sampleReplace(const std::vector<double>& prob, int n);
vint sampleNoReplace(const std::vector<double>& prob, int n);

arma::mat cov2cor(const arma::mat& sigma);
/*
double dmvnorm(const arma::colvec& x, const arma::colvec& mu, const arma::mat& sigma);
double dbvnormMarg(double k, int index, const arma::colvec& mu, const arma::mat& sigma);
double dbvnormCond(double k, double fixed, int indexFixed, const arma::colvec& mu, const arma::mat& sigma);
double dbvnormCondInterval(double k, double fixedLower, double fixedUpper, int indexFixed, const arma::colvec& mu, const arma::mat& sigma);
*/

} // namespace jmlib

#endif

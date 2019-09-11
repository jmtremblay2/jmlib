#ifndef PROBIT_UTILS
#define PROBIT_UTILS

#include <armadillo>

namespace rentix{

// turns a lower triangular matrix into a vector of its elements
arma::colvec mat2vec(const arma::mat& x);
arma::mat vec2mat(const arma::colvec& x);
arma::mat condCov(const arma::mat& sigma);
arma::colvec condMean(const arma::mat& sigma, const arma::colvec& mu, double obs);
arma::mat cov2diff(int n, int baseUt);
arma::mat reparamErrors(int n, int currentBase, int newBase);

}
#endif

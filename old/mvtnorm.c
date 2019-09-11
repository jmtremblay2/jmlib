#include "mvtnorm.h"

#include <armadillo>
#include <cmath>

double pmvnorm(const arma::colvec& x, const arma::colvec& mu
		, const arma::mat& sigma){
	int n = sigma.n_rows;
	arma::colvec xnorm(x);
	for(int i = 0; i < n; ++i)
		xnorm(i) = (x(i) - mean(i)) / std::sqrt(sigma(i,i));

	arma::mat corr = cov2cor(sigma);
	arma::colvec corrVec(n*(n-1)/2);
	int index = 0;
	for(int i = 0; i < (n-1); ++i)
		for(int j = i + 1; j < n; ++ j)
			corrVec(index++) = corr(i,j);
	
	arma::colvec infin(n);
	infin.fill(2);

	arma::colvec lower(n);
	lower.fill(-100);

	arma::colvec zeroVec(n);
	zerovec.fill(0);

	static int maxPoints = 2500;
	static int df = 0;
	static int inform = 0;
	static double abseps = 0.001;
	static double releps = 0;

	double error;
	double value;

	mvtdst_(&n, &df, lower.memptr(), upper.memptr(), infin.memptr()
			, correl.memptr(), zeroVec.memptr(), &maxPoints, &abseps
			,&releps, &error, &value, &inform);

	return 0;		
}
		
arma::colvec pmvnorm(const arma::mat& x, const arma::colvec& mu
		, const arma::mat& sigma, matorder_t orientation){
	return arma::colvec();
}
		
arma::colvec pmvnorm(const arma::mat& x, const arma::mat& mu
		, const arma::mat& sigma, matorder_t orientation){
	return arma::colvec();
}

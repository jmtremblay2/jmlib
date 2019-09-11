#include "probit_utils.hpp"
#include <cmath>

arma::colvec rentix::mat2vec(const arma::mat& x){
	if(x.n_rows != x.n_cols)
		return arma::colvec();
	int size = (x.n_rows * (x.n_rows + 1)) / 2;
	arma::colvec ret(size);
	
	int index = 0;
	for(unsigned int i = 0; i < x.n_rows; ++i)
		for(unsigned int j = 0; j <= i; ++j)
			ret(index++) = x(i,j);
	
	return ret;
	}
	
arma::mat rentix::vec2mat(const arma::colvec& x){
	int n = (int) ((-1 + std::sqrt(1 + 8 * x.n_elem)) / 2 + 0.1);
	arma::mat ret(n,n);
	
	int index = 0;
	for(int i = 0; i < n; ++i)
		for(int j = 0; j <= i; ++j)
			ret(i,j) = x(index++);
			
	return ret;
	}
	
arma::mat rentix::condCov(const arma::mat& sigma){
	int n = sigma.n_rows;
	arma::mat sigma11(n - 1, n - 1);
	for(int i = 0; i < n -1; ++i)
		for(int j = 0; j < n - 1; ++j)
			sigma11(i,j) = sigma(i,j);
	
	
	arma::colvec sigma12(n-1);
	for(int i = 0; i < n-1; ++i)
		sigma12(i) = sigma(i,n-1);

	arma::rowvec sigma21(arma::mat(sigma12).t());
	arma::mat sigma22(1,1);
	sigma22(0,0) = sigma(n-1,n-1);
	
	return sigma11 - sigma12 * sigma22.i() * sigma21;
	}
	
arma::colvec rentix::condMean(const arma::mat& sigma, const arma::colvec& mu,
		double obs){
	int n = sigma.n_rows;
	
	arma::colvec obsVec(1);
	obsVec(0) = obs;
	
	arma::colvec mu1(n-1), mu2(1);
	for(int i = 0; i < n-1; ++i)
		mu1(i) = mu(i);
	mu2(0) = mu(n-1);
	
	arma::colvec sigma12(n-1);
	for(int i = 0; i < n-1; ++i)
		sigma12(i) = sigma(i,n-1);

	arma::mat sigma22(1,1);
	sigma22(0,0) = sigma(n-1,n-1);
	
	return mu1 + sigma12 * sigma22.i() * (obsVec - mu2);
	}
	
arma::mat rentix::cov2diff(int n, int baseUt){
	int ndiff = n-1;
	arma::mat M(ndiff, n);
	M.fill(0);
	
	// fill the baseUt column with -1, to substract the utility of baseUt
	for(int i = 0; i < ndiff; ++i)
		M(i, baseUt) = -1;
	
	// put an identity matrix for the other colums
	for(int i = 0; i < ndiff; ++i)
		M(i, i + ( i >= baseUt)) = 1;

	return M;
	}
	
arma::mat rentix::reparamErrors(int n, int currentBase, int newBase){
	if(currentBase == newBase)
		return arma::eye(n-1,n-1);
	
	int minus1_col = newBase - (newBase > currentBase);
	int zero_row = currentBase - (currentBase > newBase);
	
	arma::mat M(n-1,n-1);
	M.fill(0);
	
	for(int i = 0; i < n -1; ++i)
		M(i, minus1_col) = -1;
	
	for(int i = 0; i < n-2; ++i){
		int row = i + (i >= zero_row);
		int col = i + (i >= minus1_col);
		M(row,col) = 1;
		}
		
	return M;
	}


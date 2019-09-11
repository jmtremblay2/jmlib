#include "random.hpp"
#include "constants.hpp"
#include "types.hpp"

#include <cmath>
#include <algorithm>

#include <boost/random/random_number_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/extreme_value_distribution.hpp>
#include <boost/random/uniform_01.hpp>

using namespace std;

namespace jmlib{

	double rnorm(double mu, double sigma2){
		static boost::normal_distribution<> nd(0,1); // normal dist
		static boost::random::variate_generator<boost::mt11213b&,
			boost::normal_distribution<> > var_nor(mt, nd); // output
		return std::sqrt(sigma2) * var_nor() + mu;
	}

	double dnorm(double x, double mu, double sigma2){
		static boost::math::normal normal_d;
		double sigma = std::sqrt(sigma2);
		return boost::math::pdf(normal_d, (x - mu) / sigma ) / sigma;
	}

	double pnorm(double x, double mu, double sigma2){
		static boost::math::normal normal_d;
		double sigma = std::sqrt(sigma2);
		return boost::math::cdf(normal_d, (x - mu) / sigma ) / sigma;
	}

	double rgumbel(){
		static boost::random::extreme_value_distribution<> gumbel(0,1);
		static boost::random::variate_generator<boost::mt11213b&,
			boost::random::extreme_value_distribution<> > var_gumbel(mt, gumbel);
		return var_gumbel();	
	}

	double runif(){
		static boost::random::uniform_01<> u01;
		static boost::random::variate_generator<boost::mt11213b&,
			boost::random::uniform_01<> > var_unif(mt, u01);
		return var_unif();
	}		

	vint sampleReplace(const std::vector<double>& prob, int n){
		int N = prob.size();

		vint s(n);

		double totalProb = 0;
		for(int i = 0; i < N; ++i)
			totalProb += prob.at(i);

		//cout << totalProb << endl;
		vdouble probCum(prob);
		for(int i = 1; i < N; ++i)
			probCum.at(i) += probCum.at(i-1);

		for(int i = 0; i < n; ++i){
			double u = totalProb * runif();
			//cout << "\t" << u << endl;
			//cout << u << "\n";
			if( u < probCum.at(0)){
				s.at(i) = 0;
				continue;
			}

			for(int j = 1; j < N; ++j){
				// cout << " hello  " << j << endl;
				if(probCum.at(j-1) < u && u <  probCum.at(j)){
					s.at(i) = j;
					break;
				}
			}
		}
		return s;
	}

	vint sampleNoReplace(const std::vector<double>& prob, int n){
		// inspired from R's sample functions
		int N = prob.size();

		vint s(n);

		vdouble probCum(prob);
		vint perm(N);

		for(int i = 0; i < N; ++i)
			perm.at(i) = i;
		for(int i = 1; i < N; ++i)
			probCum.at(i) += probCum.at(i-1);
		//for(int i = 0; i < N; ++i)
		//	cout << probCum.at(i) << endl;

		double totalProb = probCum.back();
		int NthisFar = N;
		for(int i = 0; i < n; ++i){
			int choice_id = -1;
			// calculate the index of the choice
			double u = totalProb * runif();

			//	cout << "sample # " << i << " u = " << u << endl;
			for(int j = 0; j < NthisFar; ++j)
				if((j == 0 && u < probCum.at(j)) 
						|| (j != 0 &&(probCum.at(j-1) < u && u < probCum.at(j)))){
					choice_id = j;
					break;
				}

			int choice = perm.at(choice_id);
			//cout << "\tchoice " << choice << endl;
			s.at(i) = choice;
			// remove the choice from the list of possible choices
			// calculate choice
			// update total prob and N
			totalProb -= prob.at(choice);
			NthisFar --;

			// discard the selected item from permutation table and cumulative prob
			// and shift the proba and perm above it 
			for(int j = choice_id; j < NthisFar; ++j){
				probCum.at(j) = probCum.at(j+1);
				probCum.at(j) -= prob.at(choice);
				perm.at(j) = perm.at(j+1);
				//cout << perm.at(j) << endl;
			}

			// perm left
			//cout << "perm left : \n";
			//for(int j = 0; j < NthisFar; ++j)
			//	cout << probCum.at(j) << "\n";
			//cout << "\n";

		}	
		return s;
	}


	arma::mat cov2cor(const arma::mat& sigma){
		arma::mat corr(sigma);
		for(unsigned int i = 0; i < sigma.n_rows; ++i){
			double cx = 1.0/std::sqrt(sigma(i,i));
			for(unsigned int j = 0; j < sigma.n_rows; ++j)
				if(j != i){
					corr(j,i) *= cx;
					corr(i,j) *= cx;
				}
			corr(i,i) = 1;	
		}

		return corr;
	}
}


/*
	 double jmlib::dmvnorm(const arma::colvec& x, const arma::colvec& mu, const arma::mat& sigma){
// check dimensions
if( ! (x.n_elem == mu.n_elem ||
mu.n_elem == sigma.n_rows ||
sigma.n_rows == sigma.n_cols))
throw "dimensions mismatch in dmvnorm";

int dim = x.n_elem;

double norm_const = std::pow(2*jmlib::pi,(-(double) dim / 2.0));
return norm_const / std::sqrt(arma::as_scalar(arma::det(sigma))) 
 * std::exp(arma::as_scalar(-0.5 * arma::trans(x - mu) * arma::inv(sigma) * (x - mu)));
 }

 double jmlib::dbvnormMarg(double k, int index, const arma::colvec& mu, const arma::mat& sigma){
// the marginal density just considers the mean and variance term separately so there is nothing to do
return jmlib::dnorm(k,mu(index),sigma(index,index));
}


double jmlib::dbvnormCond(double k, double fixed, int indexFixed, const arma::colvec& mu, const arma::mat& sigma){
double correlation = sigma(0,1) / std::sqrt(sigma(0,0) * sigma(1,1));

}

 */

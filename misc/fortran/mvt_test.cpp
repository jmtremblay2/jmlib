#include "mvtnorm.hpp"
//#include <stdio.h>
#include <iostream>
#include <time.h>
#include <jm/random.hpp>
#include <jm/mvtnorm.hpp>

using namespace std;

int main(){
	printf("hello\n");
	int n = 3;
	int df = 0;
	double lower[] = {-100,-100,-100};
	double upper[] = {-0.1443376, -0.1443376, -0.1443376};
	int infin[] = {2,2,2};
	double correl[] = {2.0/3.0, 1.0/3.0, 2.0/3.0};
	double delta[] = {0,0,0};
	int maxpts = 2500;
	double abseps = 0.001;
	double releps = 0;
	double error;
	double value;
	int inform;
/*
	mvtdst_(n,df,lower,upper,infin,correl,delta,maxpts,abseps,releps,error
			,value,inform);
*/
/*
int mvtdst_(integer *n, integer *nu, doublereal *lower, 
	doublereal *upper, integer *infin, doublereal *correl, doublereal *
	delta, integer *maxpts, doublereal *abseps, doublereal *releps, 
	doublereal *error, doublereal *value, integer *inform__);
*/

	clock_t start = clock();
/*	for(int i = 0; i < 1000; ++i)
	mvtdst_(&n,&df,lower,upper,infin,correl,delta, &maxpts, &abseps,&releps
			,&error, &value,&inform);
*/	clock_t end = clock();
	printf("%lf\n",value);
	printf("%lf\n",((double)(end-start)/CLOCKS_PER_SEC));

	arma::mat sigma("3 2 1;2 3 2;1 2 3");
	//arma::mat corr = jmlib::cov2cor(sigma);
	arma::colvec x("0.25 1.25 2.25");
	arma::colvec mean("0.5 1.5 2.5");
	double val;
	start = clock();
	for(int i = 0; i < 10000; ++i)
		val = jmlib::pmvnorm(x, mean, sigma);
	end = clock();

	cout << jmlib::pmvnorm(x, mean, sigma) << endl;
	cout << ((double)(end-start)/CLOCKS_PER_SEC) << endl;
	return 0;
}

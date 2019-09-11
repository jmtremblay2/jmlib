#include "discrete_cont1.hpp"
#include "dataut.hpp"
#include "datareg.hpp"
#include "jmlib_types.hpp"
#include "jmlib_gradient.hpp"
#include "jmlib_random.hpp"
#include "model.hpp"
#include "logit.hpp"
#include "jmlib_constants.hpp"
#include "regression.hpp"

#include <iostream>
#include <limits>

rentix::DC4::DC4(std::string datafile, std::string model) :
		Model(model),
		Ddisc(datafile,model),
		Dcont(datafile,model)
		{
	fillParams(model);	
	fillNpar();
	fillNames();
	fillBounds();
	initSolverContainers(low,up);
	genErrors(nalt);
	fillStart();	
	}


void rentix::DC4::fillParams(std::string model){	
	nalt = Ddisc.getY().max() + 1;
	nobs = Ddisc.getNobs();
	nsim = jmlib::getVal<int>("nsim",model,1000);
	return;
}

void rentix::DC4::fillNpar(){
	// parameter is : (B_disc, B_reg, Chol_line2, ... Chol_line(nalt))
	npar.push_back(Ddisc.getNvar());
	npar.push_back(Dcont.getNvar());
	for(int i = 2; i <= nalt; ++i)
		npar.push_back(i);
		
	nparTotal = 0;
	for(int i = 0; i < npar.size(); ++i)
		nparTotal += npar.at(i);
	return;
}

void rentix::DC4::fillNames(){
	// names of discrete cx
	names.push_back(jmlib::vstring(npar.at(0)));
	for(int i = 0; i < npar.at(0); ++i)
		names.at(0).at(i) = Ddisc.getVarName(i);
	// names of continuous cx
	names.push_back(jmlib::vstring(npar.at(1)));
	for(int i = 0; i < npar.at(1); ++i)
		names.at(1).at(i) = Dcont.getVarName(i);
	// the cholesky elements
	// l_0,0 = 1 (don't estimate this one)
	// l_1,0  l_1,1
	// l_2,0  l_2,1  l_2,2
	// ...
	// l_nalt-1,0 ... l_nalt-1,nalt-1 
	for(int i = 1; i < nalt; ++i){
		names.push_back(jmlib::vstring(i + 1));
		for(int j = 0; j <= i; ++j){
			std::stringstream ss;
			ss << "L_{" << i << "," << j << "}";
			names.at(i+1).at(j) = ss.str();
			}
		}
	return;
}

void rentix::DC4::fillBounds(){	
	double bound = 10;
	low = jmlib::vvdouble(npar.size());
	up = jmlib::vvdouble(npar.size());
	for(int i = 0; i < npar.size(); ++i){
		low.at(i) = jmlib::vdouble(npar.at(i), -bound);
		up.at(i) = jmlib::vdouble(npar.at(i), bound);
		}
		
	return;
}

void rentix::DC4::genErrors(int err_size) {	
	for(int obs = 0; obs < nobs; ++obs){
		errorsInd.push_back(arma::mat(err_size,nsim));
		for(int alt = 0; alt < err_size; ++alt)
			for(int sim = 0; sim < nsim; ++sim)
				errorsInd.at(obs)(alt,sim) = jmlib::rnorm();
		}
	}

void rentix::DC4::fillStart(){
	start = jmlib::vvdouble(npar.size());
	for(int i = 0; i < 2; ++i)
		start.at(i) = jmlib::vdouble(npar.at(i));
	for(int i = 2; i < npar.size(); ++i){
		start.at(i) = jmlib::vdouble(npar.at(i));
		start.at(i).back() = defaultStartSD;
		}
	return;
	}

arma::colvec rentix::DC4::probas(const jmlib::vvdouble& x) const {
	
	//get the parameters out of x
	const jmlib::vdouble& beta_disc = x.at(0);
	const jmlib::vdouble& beta_cont = x.at(1);
	arma::mat L(nalt + 1, nalt + 1);
	
	// rewrite the cholesky factor
	L.fill(0);
	L(0,0) = 1;
	for(int i = 1; i < nalt + 1; ++i)
		for(int j = 0; j <= i; ++j)
			L(i,j) = x.at(i+1).at(j);

	arma::mat V = Ddisc.getUtilities(beta_disc);
	arma::colvec Y_pred_cont = Dcont.getXBeta(beta_cont);
	arma::colvec p(nobs);
	
	// CRAP to print details of probit probabilities and regression
	// residuals at convergence
	std::ofstream out_probas;("/home/jm/probas_dc.txt");
	std::ofstream out_residuals;("/home/jm/residuals_dc.txt");
	if(this->print_res){
		out_probas.open("/home/jm/probas_dc.txt");
		out_residuals.open("/home/jm/residuals_dc.txt");		
		}
	// END CRAP	
	
	for(int obs = 0; obs < nobs; ++obs){
		double nsuccess = 0; // double for division ...
		arma::mat errorsDep = L * errorsInd.at(obs);
		int Y = Ddisc.getY(obs);
		arma::rowvec Vobs = V.row(obs);

		double mean = 0;
		double var = 0;
		for(int sim = 0; sim < nsim; ++sim){
			arma::rowvec U = Vobs + trans(errorsDep.col(sim).subvec(0,nalt-1));
			int biggest_index = 0;
			double biggest_U = U(0);
			
			for(int alt = 1; alt < nalt; ++alt)
				if(U(alt) > biggest_U){
					biggest_index = alt;
					biggest_U = U(alt);
					}
					
			if(Y == biggest_index){
				// print residuals
				if(this->print_res)
					out_residuals << errorsDep(nalt,sim) << " ";
				nsuccess++;
				mean += errorsDep(nalt,sim);
				var += errorsDep(nalt,sim) * errorsDep(nalt,sim);
				}			
			} // sim
		
		
		// if theres no succes we set num of succes to .1
		// double p_obs = -1;
		
		// CRAP to print number of success
		if(this->print_res){
			out_probas << (nsuccess / nsim) << std::endl;
			out_residuals << std::endl;				
			}
		
		double p_probit;
		double f_reg;	
		
		double p_obs;
		if(nsuccess < 0.5){
			arma::mat S = L * trans(L);
			p_obs = .1 / ((double)nsim) * jmlib::dnorm(Dcont.getY(obs),Y_pred_cont(obs),S(nalt,nalt));
			}

		else if ( nsuccess < 1.5){
			arma::mat S = L * trans(L);			
			p_obs = 1.0 / ((double)nsim) * jmlib::dnorm(Dcont.getY(obs),Y_pred_cont(obs),S(nalt,nalt));
			}
		else {
			// this is the case when we have enough success in the simulation
			mean /= nsuccess;
			var /= nsuccess;
			var -= mean*mean;
			p_obs = nsuccess / ((double) nsim) * jmlib::dnorm(Dcont.getY(obs), Y_pred_cont(obs) + mean, var); 
			}
		if(p_obs == 0)
			p_obs = std::numeric_limits<double>::epsilon();
		p(obs) = p_obs;
		
		}
	
	// CRAP close the output file
	if(this->print_res){
		out_probas.close();
		out_residuals.close();
		}
		
	return p;

	}


void rentix::DC4::solve(double def) {
	Model::solve(start);
	}

void rentix::DC4::resampleY(const jmlib::vvdouble& x){
	const jmlib::vdouble& beta_hat_probit = x.at(0);
	const jmlib::vdouble& beta_hat_reg = x.at(1);
	arma::mat L_hat(nalt+1,nalt+1);
	L_hat.fill(0);
	L_hat(0,0) = 1;
	for(int row = 1; row < (nalt+1); ++row)
		for(int col = 0; col <= row; ++col)
			L_hat(row,col) = x.at(row+1).at(col); 
	
	arma::mat U = Ddisc.getUtilities(beta_hat_probit);
	arma::colvec Y_hat = Dcont.getXBeta(beta_hat_reg);
	
	for(int obs = 0; obs < Ddisc.getNobs(); ++obs){
		arma::colvec err(nalt+1);
		for(int i = 0; i < nalt + 1; ++i)
			err(i) = jmlib::rnorm();
		arma::colvec err_corr = L_hat * err;
		
		// reset continuous Y
		//cout << Dcont.getY(obs);
		Dcont.resetY(obs, Y_hat(obs) + err_corr(nalt));
		//cout << "   " << Dcont.getY(obs) << "\n";
		
		// reset discrete Y
		arma::rowvec Uobs = U.row(obs);
		for(int i = 0; i < nalt; ++i)
			Uobs(i) += err_corr(i);
		int newY = 0;
		double biggestU = -1000000000;
		for(int i = 0; i < nalt; ++i)
			if(Uobs(i) > biggestU){
				newY = i;
				biggestU = Uobs(i);
				}
		//cout << "(" << Ddisc.getY(obs) << "," << newY << ")  ";
		Ddisc.resetY(obs,newY);
		//std::cout << newY << "  ";
		}
	//std::cout << "\n\n\n";
	}
	
jmlib::vvdouble rentix::DC4::paramH2Comp(const jmlib::vvdouble& x) const {
	jmlib::vvdouble x_reparam(x);
	arma::mat Sigma(nalt + 1, nalt + 1);
	
	Sigma.fill(0);
	Sigma(0,0) = 1;
	for(int i = 1; i < (nalt + 1); ++i)
		for(int j = 0; j <= i; ++j)
			Sigma(i,j) = x.at(i+1).at(j);
			
	arma::mat L = trans(chol(Sigma));
	for(int i = 1; i < (nalt + 1); ++i)
		for(int j = 0; j <= i; ++j)
			x_reparam.at(i+1).at(j) = L(i,j);

	return x_reparam;	
	}

jmlib::vvdouble rentix::DC1::paramComp2H(const jmlib::vvdouble& x) const {
	jmlib::vvdouble x_reparam(x);
	arma::mat L(nalt + 1,nalt + 1);
	
	// get L
	L.fill(0);
	L(0,0) = 1;
	for(int i = 1; i < (nalt + 1); ++i)
		for(int j = 0; j <= i; ++j)
			L(i,j) = x.at(i+1).at(j);
			
	// compute Sigma from L
	arma::mat Sigma = L * arma::trans(L);
	
	// write Sigma into x_reparam
	for(int i = 1; i < (nalt + 1); ++i)
		for(int j = 0; j <= i; ++j)
			x_reparam.at(i+1).at(j) = Sigma(i,j);
	return x_reparam;
	}
	
void rentix::DC1::crap(){
	double solution_hard[] = {0.0903903,	0.506901,	0.0471617,	-0.766854,	0.0310905,	
			-0.0967631,	0.368132,	-0.0902815,	0.12576,	0.111556,	0.102324,	-0.194124,
			-0.0287287,	-0.143887,	-0.187351,	-0.328078,	0.073763,	0.316336,	-0.358723,
			-0.0301723,	-0.130182,	-0.289648,	-0.593399,	0.000825114,	0.282072,	
			-0.367235,	-0.415367,	-0.264652,	1.12477,	0.0918163,	0.849966,	-0.211772,
			-0.225529,	-0.097706,	-0.234497, -0.115938,	1.46021,	0.284136,	-0.175501,	
			1.11204,	-0.034808,	-0.0661511,	-0.0271844,	0.907737,	0.0467724,	-0.0167018,	
			0.129137,	0.414147,	0.592088,	-0.709185,	-0.00209954,	0.32256,	0.466821,	
			0.0480015,	0.303429};
	int index = 0;
	for(int i = 0; i < npar.size(); ++i)
		for(int j = 0; j < npar.at(i); ++j)
			solution.at(i).at(j) = solution_hard[index++];
	std::cout << "LL = " << LL(solution) << "\n";
	//return;
	computeSE();
	
	solution = paramComp2H(solution);
	print(outputFile);

	return;
	}

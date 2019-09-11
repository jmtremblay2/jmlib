#include "model.hpp"
#include "jmlib.hpp"
#include "gradient.hpp"
#include "constants.hpp"
#include "random.hpp"

#include "solvers.hpp"

#include "stdutils.hpp"
#include "types.hpp"

#include <nlopt.hpp>
#include <armadillo>

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double rentix::modelWrapperLLNoConst(const jmlib::vdouble& x
		, jmlib::vdouble& grad, void* data){
	return modelWrapperLL(x,grad,data);
}

// this is the format the solver requires
double rentix::modelWrapperLL(const jmlib::vdouble& x, jmlib::vdouble& grad
		, const void* data){
	// data is a pointer to a model, we use it to call the LL function
	const rentix::Model* mod = static_cast<const rentix::Model*>(data);

	// convert x into a vvdouble and init a gradient container
	jmlib::vvdouble xvv(toVectorVector(x,mod->npar));
	jmlib::vvdouble gvv(mod->npar.size());
	for(unsigned int i = 0; i < mod->npar.size(); ++i)
		gvv.at(i) = jmlib::vdouble(mod->npar.at(i));

	// compute the LL
	double val;
	if(!grad.empty()){
		val = mod->LL(xvv, gvv);
		grad = std::toVector(gvv);
	}
	else
		val = mod->LL(xvv);

	// print detail of iterations if verbose is required
	//if(mod->verbose && ! grad.empty()){
	if(mod->verbose && ! grad.empty()){
		//std::cout << "eval with grad num : " << (mod->ncount_grad) 
		//	<< "     without grad :" << (mod->ncount_nograd) << "\n";
		std::cout << "x :" << x << "\n";
		std::cout << "grad :" << grad << "\n";
		std::cout << "LL : " << val << "\n";
	}

	return val;
}

double rentix::modelWrapperLL(const jmlib::vdouble& x, const void* data){
	jmlib::vdouble gradEmpty;
	return modelWrapperLL(x, gradEmpty, data);
}

double rentix::Model::LL(const jmlib::vvdouble& x, jmlib::vvdouble& grad) const {
	if( ! grad.empty()){
		jmlib::vdouble xv(std::toVector(x));
		jmlib::vdouble gv = jmlib::gradient(modelWrapperLL, xv, this
				, this->delta);
		grad = std::toVectorVector(gv, npar);
		ncount_grad++;
	}
	return LL(x);
}

double rentix::Model::LL(const jmlib::vvdouble& x) const {
	ncount_nograd++;
	arma::colvec p = probas(x);
	for(unsigned int i = 0; i < p.size(); ++i)
		if(0 == p(i))
			p(i) = 10e-20;
	// when we do non-parametric bootstrap observations
	// will be sampled 0,1,2,... times
	// so we count each observation accordingly

	//return arma::as_scalar(arma::dot(bootWeights, arma::log(p)));
	double LL = 0;
	for(int i = 0; i < nobs; ++i)
		LL += bootWeights(i) * weights(i) * std::log(p(i));
	return LL;
}


arma::colvec rentix::Model::probas(const jmlib::vvdouble& x) const {
	return arma::colvec();
}

rentix::Model::Model(std::string model){
	ncount_grad = 0;
	ncount_nograd = 0;
	logLik = 888888888888888888.888;
	solved = false;

	// read stopping criteria specified in the model file
	// (3rd argument is default values)
	rel_tol = jmlib::getVal<double>("rel_tol", model, 0.000001);
	abs_tol = jmlib::getVal<double>("abs_tol", model, 0.000001);
	max_iter = jmlib::getVal<int>("max_iter", model, 1000);
	max_time = jmlib::getVal<int>("max_time", model, 10000);

	// whether iterations should be printed
	verbose = jmlib::getVal<bool>("verbose", model, false);

	std::string se = jmlib::getVal<std::string>("se",model,"hessian");
	if(se.compare("hessian") == 0)
		SEmeth = jmlib::HESSIAN;
	else if(se.compare("bootstrap") == 0)
		SEmeth = jmlib::BOOTSTRAP;
	else if(se.compare("none") == 0)
		SEmeth = jmlib::NONE;
	else
		throw std::string("bad se argument");
	nboot = jmlib::getVal<int>("nboot",model,0);

	// figure what gradient formula to use
	delta = jmlib::getVal<double>("delta",model,1e-5);	

	outputFile = jmlib::getVal<std::string>("output", model, 
			"default_output.txt");
}

void rentix::Model::solve(const jmlib::vvdouble& start){
	// starting value
	jmlib::vdouble x = std::toVector(start);
	try{
		opt.optimize(x,logLik);
	}
	catch(std::runtime_error e){
		std::cout << "model was not solved\n";
		//throw std::string("Model was not solved");
	}
	solved = true;
	// computeSE will 
	solution = std::toVectorVector(x,npar);

	std::cout << solution << std::endl;
	computeSE();	
	print(outputFile);
	return;
}

// solve with starting value = def for all parameters
void rentix::Model::solve(double def){
	if(this->start.size() > 0)
		solve(start);
	return;

	jmlib::vvdouble startLoc;
	for(unsigned int i = 0; i < npar.size(); ++i)
		startLoc.push_back(jmlib::vdouble(npar.at(i),def));
	solve(startLoc);
	return;
}

// arguments low and up are bounds on the parameter space
void rentix::Model::initSolverContainers(const jmlib::vvdouble& low
		, const jmlib::vvdouble& up, const std::string& model){
	opt = nlopt::opt(nlopt::LD_LBFGS, nparTotal);
	opt.set_max_objective(rentix::modelWrapperLLNoConst,this);

	// set stopping criteria
	opt.set_xtol_abs(abs_tol);
	opt.set_xtol_rel(rel_tol);
	opt.set_maxtime(max_time);
	opt.set_maxeval(max_iter);

	// set bounds, if specified
	// lower
	if(low.size() == npar.size())
		opt.set_lower_bounds(std::toVector(low));
	// upper
	if(up.size() == npar.size())
		opt.set_upper_bounds(std::toVector(up));

	// initalize solution containers
	for(unsigned int i = 0; i < npar.size(); ++i){
		solution.push_back(jmlib::vdouble(npar.at(i)));
		standardError.push_back(jmlib::vdouble(npar.at(i)));
	}

	// init the boostrap weights. we need that even if we don't use
	// bootstrap estimation since 
	// LL = sum bootweight_i * log(p_i)
	// so bootweights need to be 1 to not change the LL in that case
	bootWeights = arma::vec(nobs);
	bootWeights.fill(1);

	// init the weights, the data must be initialized for that 
	// so we can't do that in the constructor	
	std::string weights = jmlib::getVal<std::string>("weights", model, "");	
	if(weights.compare("") == 0){
		this->weights = arma::vec(nobs);
		this->weights.fill(1);
	}
	else
		this->weights = D->getVar(weights);

	useObs = jmlib::vbool(nobs, true);

	// read starting values if specified
	std::string startfile = jmlib::getVal<std::string>("start",model,"hello");
	if(0 != startfile.compare("") && ! (start.size() > 0)){
		std::ifstream in(startfile.c_str());
		start = jmlib::vvdouble(npar.size());
		for(unsigned int i = 0; i < npar.size(); ++i){
			start.at(i) = jmlib::vdouble(npar.at(i));
			for(int j = 0; j < npar.at(i); ++j)
				in >> start.at(i).at(j);
		}
		in.close();
	}
}

bool rentix::Model::getSolved() const {
	return solved;
}

int rentix::Model::getNobs() const {
	return nobs;
}

int rentix::Model::getNparTotal() const {
	return nparTotal;
}

jmlib::vint rentix::Model::getNpar() const {
	return npar;
}

double rentix::Model::getMaxLL() const {
	if(!solved)
		throw("Model was not solved\n");
	return logLik;
}

jmlib::vvdouble rentix::Model::getSol() const {
	if(!solved)
		throw("Model was not solved\n");
	return solution;
}

jmlib::vvdouble rentix::Model::getSE() const {
	if(!solved)
		throw("Model was not solved\n");
	return standardError;
}

jmlib::vvstring rentix::Model::getNames() const {
	return names;
}

void rentix::Model::print(std::string outputfile){
	FILE* f = std::fopen(outputfile.c_str(),"w");
	fprintf(f, "\tvariable estimate standard_error\n");
	for(unsigned int i = 0; i < npar.size(); ++i)
		for(int j = 0; j < npar.at(i); ++j)
			fprintf(f,"\t%s &  %.4lf & %.4lf\\\\\n",names.at(i).at(j).c_str()
					, solution.at(i).at(j), standardError.at(i).at(j));

	fprintf(f, "\n");
	fprintf(f, "final LL = %lf\n", logLik);
	fprintf(f, "number of obs = %d\n", nobs);
	fclose(f);
}

// default method with hessian matrix
void rentix::Model::computeSE(){
	switch(SEmeth){
		case jmlib::HESSIAN:
			hessianSE();
			break;

		case jmlib::BOOTSTRAP:
			bootstrapSE();
			break;

		case jmlib::NONE:
			noSE();
			break;
	}
	return;
}

void rentix::Model::hessianSE(){

	jmlib::vdouble x(std::toVector(solution));

	// if -H is not invertible set the SD to -1
	arma::mat H = jmlib::hessian(modelWrapperLL, x, this , this->delta);	
	arma::colvec se(nparTotal);
	try{
		arma::mat cov = arma::inv( - H);
		se = sqrt(cov.diag());
	}
	catch(std::runtime_error e){
		std::cout << "-Hessian is not invertible\n";
		se.fill(-1);
	}

	// copy the SD into a vector<vector>
	jmlib::vdouble seVec(nparTotal);
	for(unsigned int i = 0; i < seVec.size(); ++i)
		seVec.at(i) = se(i);
	standardError = std::toVectorVector(seVec,npar);
	return;
} 


void rentix::Model::bootstrapSE(){
	this->verbose = false;
	jmlib::vvdouble solBackup = solution;	
	jmlib::vdouble start = toVector(solution);
	// var(b) = E(b^2) - E^2(b), use this formula here
	jmlib::vdouble Esquare(start.size()), squareE(start.size());
	for(unsigned int i = 0; i < Esquare.size(); ++i)
		Esquare.at(i) = squareE.at(i) = 0;

	for(int b = 0; b < nboot; ++b){
		if( ! (b%10))
			cout << "bootstrap iteration # " << b << endl;
		nonParamResample();
		jmlib::vdouble x_boot(start);

		// if for some reasons we cannot optimize the bootstrap
		// sample we just discard it and move to the next iteration
		double maxLLBoot;
		try{
			opt.optimize(x_boot, maxLLBoot);
		}
		catch(std::runtime_error e){
			b--;
			continue;
		}

		for(unsigned int i = 0; i < x_boot.size(); ++i){
			squareE.at(i) += x_boot.at(i);
			Esquare.at(i) += x_boot.at(i) * x_boot.at(i);
		}
	}
	// calculate means from sums
	for(unsigned int i = 0; i < squareE.size(); ++i){
		// Esquare just needs to be divided by nboot
		Esquare.at(i) /= nboot;
		// squareE needs to be divided by nboot to obtain E(b)
		// and then multiply by itself to obtain E^2(b)
		squareE.at(i) /= nboot;
		squareE.at(i) = squareE.at(i) * squareE.at(i);
	}

	jmlib::vdouble sd(squareE.size());
	for(unsigned int i = 0; i < sd.size(); ++i)
		sd.at(i) = sqrt(Esquare.at(i) - squareE.at(i));
	// record results
	//int index = 0;
	//for(unsigned int i = 0; i < npar.size(); ++i)
	//	for(int j = 0; j < npar.at(i); ++j)
	//		standardError.at(i).at(j) =  sd.at(index++);
	standardError = toVectorVector(sd, npar);
	solution = solBackup;
}

void rentix::Model::noSE(){
	for(unsigned int i = 0; i < npar.size(); ++i)
		for(int j = 0; j < npar.at(i); ++j)
			standardError.at(i).at(j) = -1;
}

void rentix::Model::nonParamResample(){
	bootWeights.fill(0);
	for(int i = 0; i < nobs; ++i)
		bootWeights(((int) (nobs * jmlib::runif())))++;
}

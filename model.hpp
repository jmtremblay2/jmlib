#ifndef  MODEL_CLASS
#define MODEL_CLASS

//#include "circular.hpp"

#include "data.hpp"

#include "types.hpp"

#include "jmlib.hpp"
#include "constants.hpp"
#include "gradient.hpp"

#include <string>
#include <armadillo>
#include <nlopt.hpp>

namespace rentix{


double modelWrapperLLNoConst(const jmlib::vdouble& x, jmlib::vdouble& grad
		, void* data);

double modelWrapperLL(const jmlib::vdouble& x, jmlib::vdouble& grad
		, const void* data);

// this can be casted into a func_ptr (types.hpp)
double modelWrapperLL(const jmlib::vdouble& x, const void* data);

class Model{

	friend double modelWrapperLL(const jmlib::vdouble& x, 
			jmlib::vdouble& grad, const void* data);

	public:

	/**
	 * @brief Model constructor reads general model specifications

	 * @details Constructor expects to read attribute specifiers and values
	 * in the specification files like this:\n
	 * atr1 val1\n
	 * atr2 val2\n
	 * atr3 val3\n
	 * etc\n

	 * Does not produce errors if an attribute is not recognized, just
	 * ignores it.\n
	 *
	 * list of attribututes:\n
	 * atribute (default_value) description\n
	 * rel_tol (0.000001) relative error to stop the solver\n
	 * abs_tol (0.000001) absolute error to stop the solver\n
	 * max_iter (1000) maximum number of iterations for the solver\n
	 * max_time (10000) maximum time for the solver\n
	 * verbose (0) print dedails or solver iterations (1) or not (0)\n
	 * se (hessian) method to use to compute the standard deviation (bootstrap
	 			, hessian or none)\n
	 * nboot (0)  number of bootstrap iterations\n
	 * delta (1e-5) delta (step size) to compute finite differences derivatives\n
	 * output (default_output.txt)  output file name for results\n
	 *
	 * @param model the name of the model specication file
	 */
	Model(std::string model);
	Model() {}
	virtual ~Model() {};

	/**
	 * @brief Compute log likelihood and gradient

	 * @details
	 * compute the log likelihood at x and writes the gradient at x
	 * in grad
	 *
	 * only one of (the two LL and probas) need to be implemented
	 *
	 * @param x a vvdouble of the parameters of the model
	 * @param grad a preallocated vvdouble to write the gradient
	 *
	 * @return the loglikelihood at x
	 * 
	 * @see jmlib::vvdouble
	 *
	 */
	virtual double LL(const jmlib::vvdouble& x, jmlib::vvdouble& grad) const;

	/**
	 * @brief computes only the loglikelihood
	 *
	 * @details
	 * computes the loglikelihood and lets the other functions compute
	 * the gradient (automatically) by finite differences
	 *
	 * only one of (the two LL and probas) need to be implemented
	 *
	 * @param x a vvdouble of the parameters of the model
	 * @return the loglikelihood at x
	 *
	 * @see jmlib::vvdouble
	 *
	 */
	virtual double LL(const jmlib::vvdouble& x) const;

	/**
	 * @brief computes the probability of each observation
	 *
	 * @details
	 * computes the probability of each observation so it can be used
	 * by the LL functions to compute the LL and its gradient. If you
	 * implement one of these functions, you do not need to implement
	 * probas
	 *
	 * only one of (the two LL and probas) need to be implemented
	 *
	 * @param x a vvdouble of the parameters of the model
	 * @return a colvec of size n that contain the probability of each obs.
	 *
	 */
	virtual arma::colvec probas(const jmlib::vvdouble& x) const; 

	/**
	 * @brief launches the solver with start as a starting point
	 *
	 * @param start the starting point for the solver
	 */
	virtual void solve(const jmlib::vvdouble& start);

	/**
	 * @brief lanches the solver with one repeated value as a starting point
	 *
	 * @param def the single value (default = 0)
	 */
	virtual void solve(double def = 0);

	bool getSolved() const;

	int getNobs() const;
	int getNparTotal() const;

	/**
	 * @brief get the size of each parameter of the model
	 *
	 * @details if the models depends on (p1, p2, p2) retuns
	 * the size of p1, the size of p2 and the size of p3 in a vector
	 * of ints
	 *
	 * @return the size of the parameter
	 */
	jmlib::vint getNpar() const;

	/*
	 * @brief get the estimated maximum likelihood.
	 *
	 * @details should verify with getSolved that the problem was
	 * succesfully maximized
	 *
	 * @return max LL
	 * @see rentix::Model::getSolved
	 */
	double getMaxLL() const;

	/*
	 * @brief get the estimated solution
	 *
	 * @details
	 * should verify that the model was succesfully optimized
	 *
	 * @return the estimated parameters that maximize the LL function
	 * @see rentix::Model::getSolved
	 */
	jmlib::vvdouble getSol() const;


	/*
	 * @brief get the estimated standard errors
	 *
	 * @details
	 * should verify that the model was succesfully optimized
	 *
	 * @return the estimated standard errors of the estimates
	 * @see rentix::Model::getSolved
	 */
	jmlib::vvdouble getSE() const;


	/*
	 * @brief get the parameter names
	 *
	 * @return the parameter names
	 */
	jmlib::vvstring getNames() const;

	// print to file the results
	void print(std::string outputfile);

	virtual void crap() {};

	// experimental functions
	void sampleFromPop();
	protected:

	// computeSE will call one of the fooSE methods
	virtual void computeSE();
	void hessianSE();
	void bootstrapSE();
	void noSE();

	void nonParamResample();

	Data *D;

	jmlib::vvdouble start;

	void initSolverContainers(
			const jmlib::vvdouble& low =jmlib::vvdouble(), 
			const jmlib::vvdouble& up = jmlib::vvdouble(),
			const std::string& model = "");

	double logLik;
	jmlib::vint npar;
	int nparTotal;
	bool solved;
	jmlib::vvdouble solution;
	jmlib::vvdouble standardError;
	//jmlib::vvdouble modGrad;
	jmlib::vvstring names;
	int nobs;

	// jmlib::vvdouble low, up;
	// weights for pseudo-max-LL
	//std::string weights;
	//bool useWeights;
	arma::vec weights;
	arma::vec bootWeights;

	jmlib::vbool useObs;				

	// if true, detail of iterations will be displayed
	bool verbose;

	// method to compute standard errors of the estimates
	jmlib::SEmethod SEmeth;
	int nboot;

	bool print_res;

	// stopping criteria
	double rel_tol, abs_tol;
	int max_iter, max_time;
	// file to write results in
	std::string outputFile;
	double delta; 
	// derivatives formula
	// jmlib::GradConstants gc;

	// count number of times 
	mutable int ncount_grad;
	mutable int ncount_nograd;

	nlopt::opt opt;

};

} // namespace rentix
#endif

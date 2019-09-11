#ifndef JMLIB_GRADIENT
#define JMLIB_GRADIENT

#include <armadillo>
#include "constants.hpp"

/**
 * basic utility functions
 */
namespace jmlib{

/**
 * SE evaluation methods
 */
enum SEmethod{
	HESSIAN,
	BOOTSTRAP,
	NONE
};

/**
 * calculates the gradient of a function
 *
 * @param f a pointer to a function that takes a const vdouble& 
 *		and a const void*
 * @param v the point to evaluate the gradient
 * @param data a pointer to data
 * @param d a finite difference constant
 *
 * @return the gradient of f evaluated at v
 */
vdouble gradient(func_ptr f, const vdouble& v, const void* data, double d);

/**
 * calculates the hessian matrix of a function
 *
 * @param f a pointer to a function that takes a const vdouble&
 *		and a const void*
 * @param v the point to evaluate the function
 * @param data a pointer to data
 * @param d a finite difference constant
 *
 * @return the hessian matrix of f evaluated at v
 */
arma::mat hessian(func_ptr f, const vdouble& v, const void* data, double d);

/**
 * calculate the derivative of a function at position pos
 *
 * @param f a pointer to a function
 * @param v the point to evaluate the fonction
 * @param a pointer to data
 * @param pos the index of x for the derivative
 * @param d a finite difference constant
 *
 * @return the derivative of f with respect to x_pos
 */
double deriv(func_ptr f, const vdouble& v, const void* data, int pos
		, double d);

/**
 * calculate the mixed derivative of a function
 *
 * @param f a pointer to a function
 * @param v the point to evaluate the fonction
 * @param a pointer to data
 * @param i the index of the first derivative
 * @param j the index of the second derivative
 * @param d a finite difference constant
 *
 * @return the mixed derivative of f with respect to x_i and x_j
 */
double mixedDeriv(func_ptr f, const vdouble& v, const void* data
		, int i, int j, double d);

}
#endif

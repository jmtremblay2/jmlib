#ifndef H_ARMA_UTILS
#define H_ARMA_UTILS

#include <armadillo>
#include "types.hpp"

/**
 * custom utilities for the arma namespace
 */
namespace arma{

/**
 * matrix equality
 *
 * @param a a matrix
 * @param b a matrix
 *
 * @return true if the difference between any two elems is less than delta
 */
bool operator==(const mat& a, const mat& b);

/**
 * replicate a line into a matrix
 *
 * @param r a row to be replicated
 * @param nrep the number of time to replicate r
 *
 * @return a matrix with each line equal to r
 */
mat repLines(const rowvec& r, int nrep);

/**
 * calculates the biggest index of each line of a matrix
 *
 * @param M a matrix
 *
 * @return a vector of the biggest indexes of M's lines
 */
icolvec maxIndex(const mat& M);

colvec toVector(const jmlib::vdouble& v);
jmlib::vdouble toVector(const colvec& v);
arma::imat toImat(const mat& M);
}


#endif

#ifndef CONSTANTS_JMLIB
#define CONSTANTS_JMLIB

#include "jmlib_types.hpp"

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

namespace jmlib {

/**
 * numerical approximation of sqrt(pi)/6
 */
const static double sqrtPi2over6 = 1.28254983016186;

/**
 * numerical approximation of pi
 */
const static double pi = 3.141592653589793;

}
/*
// Constants for polynomial interpolation methods for derivative estimation
struct GradConstants{

	GradConstants(const double diff[], const double cx[], int npoints){
		this->npoints = npoints;
		for(int i = 0; i < npoints; ++i){
			this->diff += diff[i];
			this->cx += cx[i];
			}
		}

	GradConstants(const GradConstants& gc, double delta): 
			diff(gc.diff), cx(gc.cx), npoints(gc.npoints), delta(delta){
		for(int i = 0; i < npoints; ++i)
			diff.at(i) *= delta;
		}

	GradConstants() {}
	
	jmlib::vdouble diff;
	jmlib::vdouble cx;
	double delta;
	int npoints;
	};

// 3 point formula
static const double diff_3points[] = {1, -1};
static const double cx_3points[] = {1.0/2.0, -1.0/2.0};
static const int npoints_3points = 2;
static const GradConstants gcThree(diff_3points, cx_3points, 	
		npoints_3points);

// 5 point formula
static const double diff_5points[] = {2,1,-1,-2};
static const double cx_5points[] = {-1.0/12.0,8.0/12.0,-8.0/12.0,1.0/12.0};
static const int npoints_5points = 4;
static const GradConstants gcFive(diff_5points, cx_5points,
		npoints_5points);

}
*/
#endif

#ifndef JMLIB_TYPES
#define JMLIB_TYPES

#include <vector>
#include <string>
#include <armadillo>
#include <map>

namespace jmlib{


typedef std::vector<double> vdouble;
typedef std::vector<vdouble> vvdouble;

typedef std::vector<int> vint;

typedef std::vector<arma::mat> vmat;

typedef std::vector<std::string> vstring;
typedef std::vector<vstring> vvstring;

typedef std::map<std::string, int> mapStringInt;

/**
 * SE evaluation methods
 */
enum SEmethod{
	HESSIAN,
	BOOTSTRAP,
	NONE
	};

typedef double (*func_ptr)(const vdouble&, const void*);
typedef double (*func_ptr_grad)(const vdouble&, vdouble&, const void*);

} // namespace jmlib

#endif

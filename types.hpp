#ifndef JMLIB_TYPES
#define JMLIB_TYPES

#include <vector>
#include <string>
#include <armadillo>
#include <map>

namespace jmlib{

// vector types
typedef std::vector<bool> vbool;

typedef std::vector<double> vdouble;
typedef std::vector<vdouble> vvdouble;

typedef std::vector<int> vint;

typedef std::vector<arma::mat> vmat;

typedef std::vector<std::string> vstring;
typedef std::vector<vstring> vvstring;

typedef std::map<std::string, int> mapStringInt;

// function pointers
typedef double (*func_ptr)(const vdouble&, const void*);
typedef vdouble (*grad_ptr)(const vdouble&, const void*);

//typedef double (*func_ptr_grad)(const vdouble&, vdouble&, const void*);

typedef double (*c_func)(double*, int, void*);
} // namespace jmlib

#endif

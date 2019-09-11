#include "arma_utils.hpp"
#include "types.hpp"

namespace arma{

bool operator==(const mat& a, const mat& b){
	// the maximum difference between two elements we may tolerate
	double tol = 0.00000001;
	if(a.n_rows != b.n_rows || a.n_cols != b.n_cols)
		return false;
	for(unsigned int i = 0; i < a.n_rows; ++i)
		for(unsigned int j = 0; j < a.n_cols; ++j)
			if(std::fabs(a(i,j) - b(i,j)) > tol)
				return false;
	return true;
}

mat repLines(const rowvec& r, int nrep){
	mat M(nrep, r.n_elem);
	for(int i = 0; i < nrep; ++i)
		M.row(i) = r;
	return M;
}

icolvec maxIndex(const mat& M){
	icolvec maxI(M.n_rows);
	for(unsigned int r = 0; r < M.n_rows; ++r){
		double max = M(r,0);
		int index = 0;

		for(unsigned int c = 1; c < M.n_cols; ++c)
			if(M(r,c) > max){
				max = M(r,c);
				index = c;
			}
		maxI(r) = index;
	}
	return maxI;
}



colvec toVector(const jmlib::vdouble& v){
	colvec ret(v.size());
	for(unsigned int i = 0; i < v.size(); ++i)
		ret(i) = v.at(i);
	return ret;
	}

jmlib::vdouble toVector(const colvec& v){
	jmlib::vdouble ret(v.n_elem);
	for(unsigned int i = 0; i < v.n_elem; ++i)
		ret.at(i) = v(i);
	return ret;
	}

arma::imat toImat(const arma::mat& M){
	arma::imat IM(M.n_rows, M.n_cols);
	for(unsigned int i = 0; i < M.n_rows; ++i)
		for(unsigned int j = 0; j < M.n_cols; ++j)
			IM(i,j) = (int) (M(i,j) + 0.5);
	return IM;
}
}

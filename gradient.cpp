#include "gradient.hpp"
#include "types.hpp"

#include <iostream>
#include <omp.h>

jmlib::vdouble jmlib::gradient(func_ptr f, const vdouble& v,  const void* data
			, double d){
	//std::cout << "in gradient func :\n";	
	vdouble grad(v.size());
//#pragma omp parallel for num_threads(12)
	for(unsigned int i = 0; i < v.size(); ++i)
		grad.at(i) = deriv(f,v,data,i,d);
	return grad;
}

arma::mat jmlib::hessian(func_ptr f, const vdouble& v, const void* data
			, double d){
	arma::mat H(v.size(), v.size());
	H.fill(0);
	for(unsigned int i = 0; i < v.size(); ++i)
		for(unsigned int j = 0; j <= i; ++j)
			H(i,j) = H(j,i) = mixedDeriv(f,v,data,i,j,d);

	return H;
}

double jmlib::deriv(func_ptr f, const vdouble& v, const void* data, int pos
			, double d){
	vdouble vmin(v), vmax(v);
	vmin.at(pos) -= d;
	vmax.at(pos) += d;
	return (f(vmax,data) - f(vmin,data)) / (2*d);
}

double jmlib::mixedDeriv(func_ptr f, const vdouble& v, const void* data
			, int i, int j, double d){
	static double cx[] = {1,-1,-1,1};
	static double diff1[] = {1,-1,1,-1};
	static double diff2[] = {1,1,-1,-1};
	double h = 0;
	for(int k = 0; k < 4; ++k){
		vdouble v2(v);
		v2.at(i) += diff1[k] * d;
		v2.at(j) += diff2[k] * d;
		h += cx[k] * f(v2,data);
		}
	double mixed = h / (4 * d * d);
	return mixed;
	}
/*
jmlib::vvdouble jmlib::gradient(const rentix::Model* m, 
		const jmlib::vvdouble& x, const jmlib::GradConstants& gc){
	
	jmlib::vvdouble grad(x.size());
	jmlib::vvdouble beta(x);
	for(int i = 0; i < x.size(); ++i){
		grad.at(i) = jmlib::vdouble(x.at(i).size());
		for(int j = 0; j < x.at(i).size(); ++j){
			double der = 0;
			for(int p = 0; p < gc.npoints; ++p){
				beta.at(i).at(j) += gc.diff.at(p);
				der += gc.cx.at(p) * m->LL(beta);
				beta.at(i).at(j) -= gc.diff.at(p);
				}
			grad.at(i).at(j) = der / gc.delta;
			}
		}
	return grad;
	}
	
// From David Eberly	
arma::mat jmlib::hessian(const rentix::Model* m, const jmlib::vvdouble& x, 
		const jmlib::GradConstants& gc){
	arma::mat H(m->nparTotal,m->nparTotal);
		
	// we compute the hessian with respect to the
	// human-readable version of the parameter
	jmlib::vdouble x_human_vec = jmlib::vv2v(m->paramComp2H(x));
	
	for(int i = 0; i < m->nparTotal; ++i)
		for(int j = 0; j <= i; ++j)
			H(i,j) = H(j,i) = jmlib::mixedDer(m,x_human_vec,gc,i,j);
	return H;
	}

// From David Eberly
double jmlib::mixedDer(const rentix::Model* m, jmlib::vdouble& x_human_vec,
		const jmlib::GradConstants& gc, int index1, int index2){
	double sum = 0;
	for(int i = 0; i < gc.npoints; ++i)
		for(int j = 0; j < gc.npoints; ++j){
			// get the finite differences
			x_human_vec.at(index1) += gc.diff.at(i);
			x_human_vec.at(index2) += gc.diff.at(j);
					
			// turn that into a machine vvdouble
			jmlib::vvdouble x_comp_vv = 
					m->paramH2Comp(jmlib::v2vv(x_human_vec, m->npar));
			
			sum += gc.cx.at(i) * gc.cx.at(j) * m->LL(x_comp_vv);
					 
			x_human_vec.at(index1) -= gc.diff.at(i);
			x_human_vec.at(index2) -= gc.diff.at(j);
			}
	
	sum /= (gc.delta * gc.delta);
	return sum;
	}

arma::mat jmlib::hessian(double (*f)(const vdouble&, vdouble&, void*),
		const vdouble& x,	void* data, const GradConstants& gc){

	int n = x.size();
	jmlib::vdouble point(x);
	jmlib::vdouble grad(n);
	arma::mat H(n,n);

	for(int i = 0; i < n; ++i){
		for(int j = 0; j <= i; ++j){
			double mixed = 0;
			for(int p = 0; p < gc.npoints; ++p){
				point.at(j) += gc.diff.at(p);
				f(point, grad, data);
				mixed += gc.cx.at(p) * grad.at(i);
				point.at(j) -= gc.diff.at(p);
				}
			H(i,j) = H(j,i) = mixed / gc.denominator;
			}
		}
	return H;
	}
*/

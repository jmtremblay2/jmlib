#ifndef JMLIB_SOLVERS
#define JMLIB_SOLVERS

#include "types.hpp"

namespace jmlib{

enum MaxOrMin{
	MAXIMIZE,
	MINIMIZE
};

class Solver{
	public:
		Solver() {}
		Solver(MaxOrMin mom
				,func_ptr f
				,grad_ptr grad
				,const void* data
				,const vdouble& start
				,double delta);

		vdouble getSol() const;
		bool optimSuccess() const;

		virtual void solve() = 0;

		void setObjective(MaxOrMin mom);
		void setFunction(func_ptr f);
		void setGrad(grad_ptr grad);
		void setData(void* data);
		void setStart(const vdouble& start);
		void setDelta(double delta);

		MaxOrMin getObjective() const { return mom; }
		double getDelta() const { return delta; }
		

		double eval(const vdouble& x) const;
	

	protected:
		MaxOrMin mom;
		func_ptr f;
		grad_ptr grad;
		const void* data;
		vdouble start;
		double delta;
		bool solved;
		int n;
		vdouble solution;
};

class BFGS: public Solver{
	public:
		BFGS(){}
		BFGS(MaxOrMin mom
				,func_ptr f
				,grad_ptr grad
				,const void* data
				,const vdouble& start
				,double delta);
		virtual void solve();
};

} // namespace
#endif

#include "solvers.hpp"
#include "types.hpp"

jmlib::Solver::Solver(MaxOrMin mom
				,func_ptr f
				,grad_ptr grad
				,const void* data
				,const vdouble& start
				,double delta):
				
				f(f)
				,grad(grad)
				,data(data)
				,start(start)
				,delta(delta)
				,solved(false)
				,n(start.size()){}

jmlib::vdouble jmlib::Solver::getSol() const{
	if(solved)
		return solution;
	return jmlib::vdouble();
}

bool jmlib::Solver::optimSuccess() const{
	return solved;
}

void jmlib::Solver::setObjective(MaxOrMin mom){
	this->mom = mom;
	solved = false;
}

void jmlib::Solver::setFunction(func_ptr f){
	this->f = f;
	solved = false;
}

void jmlib::Solver::setGrad(grad_ptr grad){
	this->f = f;
	solved = false;
}

void jmlib::Solver::setData(void* data){
	this->data = data;
	solved = false;
}

void jmlib::Solver::setStart(const vdouble& start){
	this->start = start;
	solved = false;
}

void jmlib::Solver::setDelta(double delta){
	this->delta = delta;
	solved = false;
}

double jmlib::Solver::eval(const jmlib::vdouble& x) const{
	if(NULL == this->f)
		throw("function to optimize has not been defined\n");
	return this->f(x,this->data);
}

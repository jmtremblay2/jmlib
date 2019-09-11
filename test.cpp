#include "logit.hpp"
#include "oprobit.hpp"
//#include <jm/dc2.hpp>
#include <iostream>

using namespace std;

int main(){
	// logit
	rentix::Logit L("test/Dlogit.txt","test/logit.spec");
	L.solve();
	L.print("test/resultsLogit.txt");
	
	// ordered probit
	rentix::Oprobit OP("test/Doprobit.txt","test/oprobit.spec");
	OP.solve();
	OP.print("test/testresultsOP.txt");
	
	return 0;
}

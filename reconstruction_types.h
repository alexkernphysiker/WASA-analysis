// this file is distributed under 
// MIT license
#ifndef IOAXYPGO
# define IOAXYPGO
#include <Genetic/paramfunc.h>
namespace Reconstruction{
	using namespace Genetic;
	typedef Add2<
		PolynomFunc<1,0,2>,
		Mul<Arg<0>,PolynomFunc<1,3,2>>
	> He3EnergyFRH1;
	typedef Add2<
		PolynomFunc<1,0,2>,
		Mul<Arg<0>,PolynomFunc<1,3,2>>
	> He3EnergyFRH2;
};
#endif
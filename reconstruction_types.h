// this file is distributed under 
// MIT license
#ifndef IOAXYPGO
# define IOAXYPGO
#include <Genetic/paramfunc.h>
namespace Reconstruction{
	using namespace Genetic;
	typedef Add4<
		                 PolynomFunc<1,0,3>,
		Mul<      Arg<0>,PolynomFunc<1,4,3>>,
		Mul< Sqr<Arg<0>>,PolynomFunc<1,8,3>>,
		Mul<Pow3<Arg<0>>,PolynomFunc<1,12,3>>
	> He3EnergyFRH1;
	typedef Add4<
		                 PolynomFunc<1,0,3>,
		Mul<      Arg<0>,PolynomFunc<1,4,3>>,
		Mul< Sqr<Arg<0>>,PolynomFunc<1,8,3>>,
		Mul<Pow3<Arg<0>>,PolynomFunc<1,12,3>>
	> He3EnergyFRH2;
};
#endif
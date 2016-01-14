// this file is distributed under 
// MIT license
#ifndef IOAXYPGO
# define IOAXYPGO
#include <Genetic/paramfunc.h>
namespace Reconstruction{
	using namespace Genetic;
	typedef Add3<
		Add3<Par<0>,Mul<Par<1>,Arg<1>>,Mul<Par<2>,Sqr<Arg<1>>>>,
		Mul< Arg<0>, Add3<Par<3>,Mul<Par<4>,Arg<1>>,Mul<Par<5>,Sqr<Arg<1>>>>>,
		Mul< Sqr<Arg<0>>, Add3<Par<6>,Mul<Par<7>,Arg<1>>,Mul<Par<8>,Sqr<Arg<1>>>>>
	> He3EnergyFRH1;
	typedef Add3<
		Add3<Par<0>,Mul<Par<1>,Arg<1>>,Mul<Par<2>,Sqr<Arg<1>>>>,
		Mul< Arg<0>, Add3<Par<3>,Mul<Par<4>,Arg<1>>,Mul<Par<5>,Sqr<Arg<1>>>>>,
		Mul< Sqr<Arg<0>>, Add3<Par<6>,Mul<Par<7>,Arg<1>>,Mul<Par<8>,Sqr<Arg<1>>>>>
	> He3EnergyFRH2;
};
#endif
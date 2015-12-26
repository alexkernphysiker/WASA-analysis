// this file is distributed under 
// MIT license
#ifndef IOAXYPGO
# define IOAXYPGO
#include <Genetic/paramfunc.h>
namespace Reconstruction{
	using namespace Genetic;
	//Lovely ROOT headers :)
	typedef Genetic::Add<
		Mul<Arg<0>,Genetic::Add<Mul<Arg<1>,Par<0>>,Par<1>>>,
		Genetic::Add<Mul<Arg<1>,Par<2>>,Par<3>>
	> He3EnergyFRH1;
	typedef Genetic::Add<
		Mul<Arg<0>,Genetic::Add<Mul<Arg<1>,Par<0>>,Par<1>>>,
		Genetic::Add<Mul<Arg<1>,Par<2>>,Par<3>>
	> He3EnergyFRH2;
};
#endif
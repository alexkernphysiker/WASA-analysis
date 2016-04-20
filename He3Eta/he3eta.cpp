// this file is distributed under 
// MIT license
#include "he3eta.h"
using namespace MathTemplates;
const Reaction&he3eta(){
	static Reaction main_react(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	return main_react;
}
const LinearInterpolation<double>&he3eta_sigma(){
	static LinearInterpolation<double> res{
		//http://arxiv.org/pdf/nucl-ex/0701072v1
		point<double>(-0.5,0.0),
		point<double>(0.0,100.0),
		point<double>(0.5,380.0),
		point<double>(1.5,400.0),
		point<double>(6.0,390.0),
		point<double>(12.0,380.0),
		//Extrapolating
		point<double>(24.0,360.0),
		point<double>(36.0,340.0)
	};
	return res;
}
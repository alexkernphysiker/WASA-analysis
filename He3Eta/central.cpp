// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Genetic/fit.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Kinematics/particles.h>
#include <Kinematics/reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
const Reaction&main_reaction(){
	static Reaction main_react(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	return main_react;
}
int main(){
	LinearInterpolation<double> sigmaHe3eta{
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
}
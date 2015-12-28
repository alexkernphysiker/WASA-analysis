// this file is distributed under 
// MIT license
#include <unistd.h>
#include <Genetic/fit.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
#include <reconstruction_fit.h>
#include "reconstruction_types.h"
using namespace std;
using namespace Genetic;
using namespace Reconstruction;
RANDOM engine;
int main(int,char**){
	Plotter::Instance().SetOutput(SimulationDataPath(),"he3reconstruction");
	ProcessFit<He3EnergyFRH1>(
		"He3.E.FRH1",
		make_shared<GenerateByGauss>()<<make_pair(0,1)<<make_pair(0,1)<<make_pair(0,1)<<make_pair(0,1),
		make_shared<Filter>([](const ParamSet&){return true;}),
		engine
	);
	
}
// this file is distributed under 
// MIT license
#include <unistd.h>
#include <Genetic/fit.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
#include <reconstruction_fit.h>
#include <str_get.h>
#include "reconstruction_types.h"
using namespace std;
using namespace Genetic;
using namespace Reconstruction;
RANDOM engine;
int main(int,char**){
  	Plotter::Instance().SetOutput(SimulationDataProcess::SimulationDataPath(),"He3Ekin");
	Plotter::Instance()<<"set xrange [0:0.4]"<<"set yrange [0:0.6]";
	SimulationDataProcess::ProcessEnergyThetaFit<He3EnergyFRH1>("He3.E.FRH1",make_pair(0,0.6),
		make_shared<GenerateByGauss>()
			<<make_pair(0,100)<<make_pair(0,100)<<make_pair(0,100)
			<<make_pair(0,100)<<make_pair(0,100)<<make_pair(0,100)
			<<make_pair(0,100)<<make_pair(0,100)<<make_pair(0,100)
		,make_shared<Filter>([](const ParamSet&){return true;}),
		engine
	);
	Plotter::Instance()<<"set xrange [0.2:0.6]"<<"set yrange [0.2:0.6]";
	SimulationDataProcess::ProcessEnergyThetaFit<He3EnergyFRH2>("He3.E.FRH2",make_pair(0.25,0.6),
		make_shared<GenerateByGauss>()
			<<make_pair(0,100)<<make_pair(0,100)<<make_pair(0,100)
			<<make_pair(0,100)<<make_pair(0,100)<<make_pair(0,100)
			<<make_pair(0,100)<<make_pair(0,100)<<make_pair(0,100)
		,make_shared<Filter>([](const ParamSet&){return true;}),
		engine
	);
	
}
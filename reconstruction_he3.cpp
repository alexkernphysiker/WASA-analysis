// this file is distributed under 
// MIT license
#include <unistd.h>
#include <gnuplot_wrap.h>
#include <Genetic/fit.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
#include <reconstruction_fit.h>
#include <str_get.h>
#include <reconstruction_types.h>
using namespace std;
using namespace Genetic;
using namespace Reconstruction;
using namespace GnuplotWrap;
RANDOM engine;
int main(int,char**){
  	Plotter::Instance().SetOutput(SimulationDataProcess::SimulationDataPath(),"He3Ekin");
	SimulationDataProcess::He3ForEtaFit<He3EnergyFRH1>("He3.E.FRH1",make_pair(0,0.6),
		make_shared<GenerateByGauss>()
			<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.01)
			<<make_pair(1,0.5)<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.01)
			<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.01)
			<<make_pair(0,0.01)<<make_pair(0,0.01)<<make_pair(0,0.01)<<make_pair(0,0.01)
		,engine
	);
	SimulationDataProcess::He3ForEtaFit<He3EnergyFRH2>("He3.E.FRH2",make_pair(0.3,0.6),
		make_shared<GenerateByGauss>()
			<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.01)
			<<make_pair(1,0.5)<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.01)
			<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.01)
			<<make_pair(0,0.01)<<make_pair(0,0.01)<<make_pair(0,0.01)<<make_pair(0,0.01)
		,engine
	);
	
}
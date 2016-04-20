// this file is distributed under 
// MIT license
#include <unistd.h>
#include <gnuplot_wrap.h>
#include <math_h/hists.h>
#include <Genetic/initialconditions.h>
#include <ReconstructionFit/reconstruction_types.h>
#include <ReconstructionFit/reconstruction_fit.h>
#include <Experiment/str_get.h>
using namespace std;
using namespace Genetic;
using namespace Reconstruction;
using namespace MathTemplates;
using namespace GnuplotWrap;
RANDOM engine;
int main(){
  	Plotter::Instance().SetOutput(SimulationDataProcess::SimulationDataPath(),"He3Ekin");
	SimulationDataProcess::ForwardEkinReconstructionFit(
		"He3.E.FRH1",make_shared<He3EnergyFRH1>(),
		BinsByStep(0.0,0.005,0.3),
		BinsByStep(0.2,0.005,0.5),
		make_shared<GenerateByGauss>()
			<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.1)
			<<make_pair(1,0.5)<<make_pair(0,0.1)<<make_pair(0,0.1)
		,engine
	);
	SimulationDataProcess::ForwardEkinReconstructionFit(
		"He3.E.FRH2",make_shared<He3EnergyFRH2>(),
		BinsByStep(0.25,0.005,0.45),
		BinsByStep(0.45,0.005,0.60),
		make_shared<GenerateByGauss>()
			<<make_pair(0,0.1)<<make_pair(0,0.1)<<make_pair(0,0.1)
			<<make_pair(1,0.5)<<make_pair(0,0.1)<<make_pair(0,0.1)
		,engine
	);
	
}
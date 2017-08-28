// this file is distributed under
// GPL license
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
int main()
{
    Plotter<>::Instance().SetOutput(ENV(OUTPUT_PLOTS),"He3-reconstruction");

    SimulationDataProcess::ForwardEkinReconstructionFit(
        "He3.E.FRH1", make_shared<He3EnergyFRH1>(),
        BinsByStep(0.0, 0.005, 0.3),
        BinsByStep(0.2, 0.005, 0.5),
        make_shared<InitialDistributions>()
        << make_shared<DistribGauss>(0, 0.1) << make_shared<DistribGauss>(0, 0.1)
        << make_shared<DistribGauss>(0, 0.1) << make_shared<DistribGauss>(1, 0.5)
        << make_shared<DistribGauss>(0, 0.1) << make_shared<DistribGauss>(0, 0.1)
        , engine
    );
}

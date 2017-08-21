// this file is distributed under 
// GPL license
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
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_central_6gamma");
    vector<string> histpath_central_reconstr={"Histograms","CentralGammas6"};
    vector<string> reaction={"bound1-6g","He3eta","He3pi0","He3pi0pi0","He3pi0pi0pi0"};
    const auto runs=PresentRuns("C");
    const string runmsg=to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs";


}


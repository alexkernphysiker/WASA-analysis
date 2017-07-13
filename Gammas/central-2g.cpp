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
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_central_2gamma");
    vector<string> histpath_central_reconstr={"Histograms","CentralGammas"};
    vector<string> reaction={"He3eta","He3pi0","He3pi0pi0","He3pi0pi0pi0"};
    const auto runs=PresentRuns("R");
    const string runmsg=to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs";

    Plot<double>()
    .Line(Hist(MC,"He3eta",histpath_central_reconstr,"GammaCount").toLine(),"3He+eta")
    .Line(Hist(MC,"He3pi0",histpath_central_reconstr,"GammaCount").toLine(),"3He+pi0")
    .Line(Hist(MC,"He3pi0pi0",histpath_central_reconstr,"GammaCount").toLine(),"3He+2pi0")
    .Line(Hist(MC,"He3pi0pi0pi0",histpath_central_reconstr,"GammaCount").toLine(),"3He+3pi0")
    <<"set title 'Gamma quanta count, Monte Carlo'"<<"set key on";
    Plot<double>()
    .Hist(Hist(DATA,"C",histpath_central_reconstr,"GammaCount"))
    <<"set title 'Gamma quanta count, DATA ("+runmsg+")'";

    Plot<double>()
    .Hist(Hist(MC,"He3eta",histpath_central_reconstr,"MMass2GammaBefore-AllBins"),"3He+eta")
    .Hist(Hist(MC,"He3pi0",histpath_central_reconstr,"MMass2GammaBefore-AllBins"),"3He+pi0")
    .Hist(Hist(MC,"He3pi0pi0",histpath_central_reconstr,"MMass2GammaBefore-AllBins"),"3He+2pi0")
    .Hist(Hist(MC,"He3pi0pi0pi0",histpath_central_reconstr,"MMass2GammaBefore-AllBins"),"3He+3pi0")
    <<"set title 'Gamma pair missing mass, Monte Carlo'"<<"set key on"<<"set xrange [2:3]";
    Plot<double>()
    .Hist(Hist(DATA,"C",histpath_central_reconstr,"MMass2GammaBefore-AllBins"))
    <<"set title 'Gamma pair missing mass, DATA ("+runmsg+")'"<<"set xrange [2:3]";
    Plot<double>()
    .Hist(Hist(MC,"He3eta",histpath_central_reconstr,"InvMass2GammaBefore-AllBins"),"3He+eta")
    .Hist(Hist(MC,"He3pi0",histpath_central_reconstr,"InvMass2GammaBefore-AllBins"),"3He+pi0")
    .Hist(Hist(MC,"He3pi0pi0",histpath_central_reconstr,"InvMass2GammaBefore-AllBins"),"3He+2pi0")
    .Hist(Hist(MC,"He3pi0pi0pi0",histpath_central_reconstr,"InvMass2GammaBefore-AllBins"),"3He+3pi0")
    <<"set title 'Gamma pair invariant mass, Monte Carlo'"<<"set key on"<<"set xrange [0.0:0.8]";
    Plot<double>()
    .Hist(Hist(DATA,"C",histpath_central_reconstr,"InvMass2GammaBefore-AllBins"))
    <<"set title 'Gamma pair invariant mass, DATA ("+runmsg+")'"<<"set xrange [0.0:0.8]";


    vector<hist<double>> norm;
    Plot<double> norm_plot;
    for(const string& r:reaction){
	const auto H=Hist(MC,r,histpath_central_reconstr,"0-Reference");
	norm.push_back(H);
	norm_plot.Hist(H,r);
    }
    norm_plot<<"set key on";
}

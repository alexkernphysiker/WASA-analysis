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
#include "he3eta.h"
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_central_6gamma");
    vector<string> histpath_central_reconstr={"Histograms","CentralGammas"};
    vector<string> reaction={"He3eta","He3pi0pi0pi0"};
    const auto runs=PresentRuns("");
    const string runmsg=to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs";

    Plot<double>()
    .Hist(Hist(MC,"He3eta",histpath_central_reconstr,"GammaEnergy"),"3He+eta")
    .Hist(Hist(MC,"He3pi0pi0pi0",histpath_central_reconstr,"GammaEnergy"),"3He+3pi0")
    <<"set title 'Gamma quanta energy, Monte Carlo'"<<"set key on";
    Plot<double>()
    .Hist(Hist(DATA,"",histpath_central_reconstr,"GammaEnergy"))
    <<"set title 'Gamma quanta energy, DATA ("+runmsg+")'";

    Plot<double>()
    .Line(Hist(MC,"He3eta",histpath_central_reconstr,"GammaCount").toLine(),"3He+eta")
    .Line(Hist(MC,"He3pi0pi0pi0",histpath_central_reconstr,"GammaCount").toLine(),"3He+3pi0")
    <<"set title 'Gamma quanta count, Monte Carlo'"<<"set key on";
    Plot<double>()
    .Hist(Hist(DATA,"",histpath_central_reconstr,"GammaCount"))
    <<"set title 'Gamma quanta count, DATA ("+runmsg+")'";


    Plot<double>()
    .Hist(Hist(MC,"He3eta",histpath_central_reconstr,"GammaEnergy6After"),"3He+eta")
    .Hist(Hist(MC,"He3pi0pi0pi0",histpath_central_reconstr,"GammaEnergy6After"),"3He+3pi0")
    <<"set title '6 Gamma total energy, Monte Carlo'"<<"set key on";
    Plot<double>()
    .Hist(Hist(DATA,"",histpath_central_reconstr,"GammaEnergy6After"))
    <<"set title '6 Gamma total energy, DATA ("+runmsg+")'";
    Plot<double>()
    .Hist(Hist(MC,"He3eta",histpath_central_reconstr,"InvMass3PairsAfter-AllBins"),"3He+eta")
    .Hist(Hist(MC,"He3pi0pi0pi0",histpath_central_reconstr,"InvMass3PairsAfter-AllBins"),"3He+3pi0")
    <<"set title '6 Gamma invariant mass, Monte Carlo'"<<"set key on"<<"set xrange [0.4:0.7]";
    Plot<double>()
    .Hist(Hist(DATA,"",histpath_central_reconstr,"InvMass3PairsAfter-AllBins"))
    <<"set title '6 Gamma invariant mass, DATA ("+runmsg+")'"<<"set xrange [0.4:0.7]";

    
    vector<hist<double>> norm;
    {
	Plot<double> mc_ncd,ref;
	mc_ncd<<"set key on";
	ref<<"set key on";
	for(const string& r:reaction){
	    const auto N=Hist(MC,r,histpath_central_reconstr,"0-Reference");
	    mc_ncd.Hist(Hist(MC,r,histpath_central_reconstr,"GammaCount")/N.TotalSum(),r);
	    ref.Hist(N,r);
	    norm.push_back(N);
	}
	mc_ncd.Hist(
	    Hist(DATA,"",histpath_central_reconstr,"GammaCount")
	    /
	    Hist(DATA,"",histpath_central_reconstr,"0-Reference").TotalSum()
	,"DATA");
    }
    for(size_t bin_num=0,bin_count=norm[0].size();bin_num<bin_count;bin_num++){
	auto Q=norm[0][bin_num].X();
	string Qmsg="Q in ["+to_string(norm[0][bin_num].X().min())+":"+to_string(norm[0][bin_num].X().max())+"] MeV";
	Plot<double>()
	.Hist(Hist(DATA,"",histpath_central_reconstr,string("InvMass3PairsAfter-Bin-")+to_string(bin_num)),"DATA")
	<<"set key on"<<"set xlabel 'MM, GeV'"<<"set xrange [0.4:0.7]"<< "set yrange [0:]"
	<< "set title '3pi0 inv. mass ("+runmsg+", "+Qmsg+")'";
	Plot<double>()
	.Hist(Hist(MC,reaction[0],histpath_central_reconstr,string("InvMass3PairsAfter-Bin-")+to_string(bin_num))/norm[0][bin_num].Y(),reaction[0])
	.Hist(Hist(MC,reaction[1],histpath_central_reconstr,string("InvMass3PairsAfter-Bin-")+to_string(bin_num))/norm[1][bin_num].Y(),reaction[1])
	<<"set key on"<<"set xlabel 'MM, GeV'"<<"set xrange [0.4:0.7]"<< "set yrange [0:]"
	<< "set title '3pi0 inv. mass (Monte Carlo, "+Qmsg+")'";
    }
}


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
    vector<string> reaction={"He3eta6g","He3pi06g"};
    vector<hist<double>> norm;
    Plotter::Instance()<<"set log y";
    {
	Plot<double> mc_ncd;
	for(const string& r:reaction){
	    mc_ncd.Hist(Hist(MC,r,histpath_central_reconstr,"neutral_tracks_count"));
	    norm.push_back(Hist(MC,r,histpath_central_reconstr,"0-Reference"));
	}
	mc_ncd<<"set key on";
	Plot<double>().Hist(Hist(DATA,"",histpath_central_reconstr,"neutral_tracks_count"))<<"set key on";
    }
    Plotter::Instance()<<"unset log y";
    Plot<double>()
    .Hist(norm[0],"3He+eta")
    .Hist(norm[1],"3He+3pi0")
    <<"set key on";
    for(size_t bin_num=0,bin_count=norm[0].size();bin_num<bin_count;bin_num++){
	auto Q=norm[0][bin_num].X();
	string Qmsg="Q in ["+to_string(norm[0][bin_num].X().min())+":"+to_string(norm[0][bin_num].X().max())+"] MeV";
	Plot<double>()
	.Hist(Hist(DATA,"",histpath_central_reconstr,string("InvMass6Gamma-Bin-")+to_string(bin_num)),"DATA")
	<<"set log y"<<"set key on"<<"set title '"+Qmsg+"'"<<"set xlabel 'MM, GeV'";
	Plot<double>()
	.Hist(Hist(MC,"He3eta",histpath_central_reconstr,string("InvMass6Gamma-Bin-")+to_string(bin_num))/norm[0][bin_num].Y(),"3He+eta")
	.Hist(Hist(MC,"He3pi0",histpath_central_reconstr,string("InvMass6Gamma-Bin-")+to_string(bin_num))/norm[1][bin_num].Y(),"3He+3pi0")
	<<"set log y"<<"set key on"<<"set title '"+Qmsg+"'"<<"set xlabel 'MM, GeV'";
    }
}


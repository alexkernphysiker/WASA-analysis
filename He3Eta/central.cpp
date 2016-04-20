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
#include "he3eta.h"
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_central");
	vector<string> histpath_central_reconstr={"Histograms","CentralGammas"};
	vector<string> reaction={"He3eta","He3pi0","He3pi0pi0","He3pi0pi0pi0"};
	vector<hist<double>> norm;
	for(const string& r:reaction)norm.push_back(Hist(MC,r,histpath_central_reconstr,"0-Reference"));
	Plot<double>().Hist(norm[0],"Simulated events");
	for(size_t bin_num=0,bin_count=norm[0].size();bin_num<bin_count;bin_num++){
		auto Q=norm[0][bin_num].X();
		string Qmsg="Q in ["+to_string(norm[0][bin_num].X().min())+":"+to_string(norm[0][bin_num].X().max())+"] MeV";
		
		hist<double> data=Hist(DATA,"",histpath_central_reconstr,string("InvMass2Gamma-Bin-")+to_string(bin_num));
		Plot<double>().Hist(data);
	}
}
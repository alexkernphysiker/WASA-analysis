// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <math_h/functions.h>
#include <math_h/error.h>
#include <Genetic/fit.h>
#include <Genetic/paramfunc.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
#include <str_get.h>
#include <gethist.h>
#include <theory.h>
#include "phys_constants.h"
#include "kinematics.h"
using namespace std;
using namespace ROOT_data;
using namespace GnuplotWrap;
int main(int,char**){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
	vector<string> histpath_forward={"Histograms","He3Forward_Reconstruction"};
	hist mc_norm_fg(MC,"He3eta",histpath_forward,"0-Reference");
	hist mc_norm_bg1(MC,"He3pi0pi0",histpath_forward,"0-Reference");
	hist mc_norm_bg2(MC,"He3pi0pi0pi0",histpath_forward,"0-Reference");
	{
		double MC_events_count=0;for(auto&p: mc_norm_fg)MC_events_count+=p.y;
		cout<<"Montecarlo evens count of "<<MC_events_count<<" detected."<<endl;
	}
	hist mc_accepted_fg(MC,"He3eta",histpath_forward,"5-Kinematic cut");
	if(mc_accepted_fg.count()!=mc_norm_fg.count())
		throw Exception<decltype(main)>("mc hists size mismatch");
	{
		hist n1(MC,"He3eta",histpath_forward,"1-AllTracks");
		hist n2(MC,"He3eta",histpath_forward,"2-FPC");
		hist n3(MC,"He3eta",histpath_forward,"3-AllCuts");
		hist n4(MC,"He3eta",histpath_forward,"4-Reconstructed");
		PlotHist().Hist("All MC events",mc_norm_fg).Hist("Forward Tracks",n1).Hist("FPC",n2).Hist("E_{dep} cuts",n3).Hist("Kinematic cuts",mc_accepted_fg)
			<<"set yrange [0:]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	}
	for(size_t bin_num=0;bin_num<mc_norm_fg.count();bin_num++)if(mc_accepted_fg[bin_num].y>1){
		hist foreground(MC,"He3eta",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist background1(MC,"He3pi0pi0",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist background2(MC,"He3pi0pi0pi0",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist measured(DATA,"He3",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		foreground/=mc_norm_fg[bin_num].Y();
		background1/=mc_norm_bg1[bin_num].Y();
		background2/=mc_norm_bg2[bin_num].Y();
		double measured_events=0;for(auto&p:measured)measured_events+=p.y;
		//ToDo:make fit
	}
}
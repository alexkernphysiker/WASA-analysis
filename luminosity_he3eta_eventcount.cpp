// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <math_h/functions.h>
#include <math_h/error.h>
#include <Genetic/fit.h>
#include <Genetic/equation.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
#include <str_get.h>
#include <gethist.h>
#include <particles.h>
#include <reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace GnuplotWrap;
using namespace Genetic;
int main(int,char**){
	RANDOM engine;
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
	vector<string> histpath_forward={"Histograms","He3Forward_Reconstruction"};
	hist mc_norm_fg(MC,"He3eta",histpath_forward,"0-Reference");
	hist mc_norm_bg1(MC,"He3pi0pi0",histpath_forward,"0-Reference");
	hist mc_norm_bg2(MC,"He3pi0pi0pi0",histpath_forward,"0-Reference");
	{// debug messaging
		cout<<"Montecarlo evens count of "<<mc_norm_fg.Total()<<" detected."<<endl;
		hist n1(MC,"He3eta",histpath_forward,"1-AllTracks");
		hist n2(MC,"He3eta",histpath_forward,"2-FPC");
		hist n3(MC,"He3eta",histpath_forward,"3-AllCuts");
		hist n4(MC,"He3eta",histpath_forward,"4-Reconstructed");
		hist n5(MC,"He3eta",histpath_forward,"5-Kinematic cut");
		PlotHist().Hist(mc_norm_fg,"All events").Hist(n1,"Forward Tracks")
		.Hist(n2,"Signal in FPC").Hist(n3,"E_{dep} cuts").Hist(n5,"Kinematic cuts")
		<<"set yrange [0:]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	}
	auto luminosity=mc_norm_fg.CloneEmptyBins()
		,acceptance_fg=mc_norm_fg.CloneEmptyBins();
	for(size_t bin_num=0;bin_num<mc_norm_fg.count();bin_num++){
		hist foreground(MC,"He3eta",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist background1(MC,"He3pi0pi0",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist background2(MC,"He3pi0pi0pi0",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist measured(DATA,"He3",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		value<double> fg_count(foreground.Total());
		if(fg_count.val()>1){
			acceptance_fg[bin_num].varY()=fg_count/mc_norm_fg[bin_num].Y();
			value<double> bg1_count(background1.Total()),bg2_count(background2.Total());
			foreground/=fg_count;
			background1/=bg1_count;
			background2/=bg2_count;
			hist model_spectrum=
				 (foreground*sigmaHe3eta(mc_norm_fg[bin_num].X().val()))
				+(background1*sigmaHe3pi0pi0(mc_norm_fg[bin_num].X().val()))
				+(background2*sigmaHe3pi0pi0pi0(mc_norm_fg[bin_num].X().val()));
			PlotHist().Hist(model_spectrum,"theory Q="+to_string(mc_norm_fg[bin_num].X().val()));
			PlotHist().Hist(measured,"data Q="+to_string(mc_norm_fg[bin_num].X().val()));
			//ToDo: calculate error for model spectrum
			value<double> experimental_count(measured.Total()),model_count(model_spectrum.Total());
			luminosity[bin_num].varY()=experimental_count/model_count;
		}else{
			acceptance_fg[bin_num].varY()=value<double>(0,0);
			luminosity[bin_num].varY()=value<double>(0,0);
		}
	}
	PlotHist().Hist(acceptance_fg,"^3He+eta")<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
	PlotHist().Hist(luminosity)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integral luminosity, nb^{-1}'";
}
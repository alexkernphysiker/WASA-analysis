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
#include <he3.h>
#include "phys_constants.h"
using namespace std;
using namespace ROOT_data;
using namespace GnuplotWrap;
using namespace Genetic;
using namespace Theory;
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
		PlotHist().Hist("All MC events",mc_norm_fg).Hist("Forward Tracks",n1).Hist("FPC",n2).Hist("E_{dep} cuts",n3).Hist("Kinematic cuts",n5)
			<<"set yrange [0:]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	}
	auto events_fg=mc_norm_fg.CloneEmptyBins(),acceptance_fg=mc_norm_fg.CloneEmptyBins();
	for(size_t bin_num=0;bin_num<mc_norm_fg.count();bin_num++){
		hist foreground(MC,"He3eta",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist background1(MC,"He3pi0pi0",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist background2(MC,"He3pi0pi0pi0",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist measured(DATA,"He3",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		value norm(foreground.Total());
		if(norm.val()>1){
			foreground/=norm;			
			acceptance_fg[bin_num].varY()=norm/mc_norm_fg[bin_num].Y();
			background1/=value(background1.Total());
			background2/=value(background2.Total());
			
			Equation<DifferentialMutations<Parabolic>> fit([&measured,&foreground,&background1,&background2](const ParamSet&P)->double{
				return ChiSq(measured,(foreground*value(P[0],0))+(background1*value(P[1],0))+(background2*value(P[2],0)),3);
			});
			fit.SetFilter([](const ParamSet&P)->bool{return (P[0]>=0)&&(P[1]>=0)&&(P[2]>=0);});
			double total_events=measured.Total();
			fit.Init(100,make_shared<GenerateByGauss>()<<make_pair(total_events,total_events)<<make_pair(0,total_events)<<make_pair(0,total_events),engine);
			cout<<"Population:"<<fit.PopulationSize()<<endl;
			cout<<"Parameters:"<<fit.ParamCount()<<endl;
			while(!fit.AbsoluteOptimalityExitCondition(0.0000001))
				fit.Iterate(engine);
			cout<<fit.iteration_count()<<" iterations"<<endl;
			cout<<"Chi^2 = "<<fit.Optimality()<<endl;
			cout<<"Parameters:"<<endl;
			cout<<fit<<endl;
			cout<<"Errors:"<<endl;
			auto errors=fit.GetParamParabolicErrors({1,1,1});
			cout<<errors<<endl;
			PlotHist().Hist(string("Data-")+to_string(bin_num),measured).Hist(string("He3eta-")+to_string(bin_num),foreground*fit[0])
			.Hist(string("He3pi0pi0-")+to_string(bin_num),background1*fit[1]).Hist(string("He3pi0pi0pi0-")+to_string(bin_num),background2*fit[2])
			<<"set xlabel 'Missing mass, GeV'"<<"set ylabel 'Events count'";
			events_fg[bin_num].varY()=value(fit[0],errors[0]);
		}
	}
	PlotHist().Hist("Acceptance",acceptance_fg)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
	PlotHist().Hist("Events He3eta",events_fg)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	auto luminosity=events_fg/(acceptance_fg*sigmaHe3eta);
	PlotHist().Hist("Luminosity He3eta",events_fg)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integral luminosity, nb^{-1}'";
}
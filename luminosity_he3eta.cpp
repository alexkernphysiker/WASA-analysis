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
#include <theory.h>
#include "phys_constants.h"
#include "kinematics.h"
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
		double MC_events_count=0;for(auto&p: mc_norm_fg)MC_events_count+=p.Y().val();
		cout<<"Montecarlo evens count of "<<MC_events_count<<" detected."<<endl;
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
		double norm=0;for(auto&p:foreground)norm+=p.Y().val();
		if(norm>0){
			foreground/=value(norm);
			
			acceptance_fg[bin_num].varY()=value(norm)/mc_norm_fg[bin_num].Y();
			double norm_b1=0;for(auto&p:background1)norm_b1+=p.Y().val();
			background1/=value(norm_b1);
			double norm_b2=0;for(auto&p:background2)norm_b2+=p.Y().val();
			background2/=value(norm_b2);
			
			double total_events=0;for(auto&p:measured)total_events+=p.Y().val();
			
			Equation<DifferentialMutations<Parabolic>> fit([&measured,&foreground,&background1,&background2](const ParamSet&P)->double{
				return ChiSq(measured,(foreground*value(P[0],0))+(background1*value(P[1],0))+(background2*value(P[2],0)),3);
			});
			fit.SetFilter([](const ParamSet&P)->bool{return (P[0]>=0)&&(P[1]>=0)&&(P[2]>=0);});
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
			.Hist(string("He3pi0pi0-")+to_string(bin_num),background1*fit[1]).Hist(string("He3pi0pi0pi0-")+to_string(bin_num),background2*fit[2]);
			events_fg[bin_num].varY()=value(fit[0],errors[0]);
		}
	}
	PlotHist().Hist("Acceptance",acceptance_fg);
	PlotHist().Hist("Events He3eta",events_fg);
}
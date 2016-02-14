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
	auto events_fg=mc_norm_fg.CloneEmptyBins(),acceptance_fg=mc_norm_fg.CloneEmptyBins(),chi_sq=mc_norm_fg.CloneEmptyBins();
	for(size_t bin_num=0;bin_num<mc_norm_fg.count();bin_num++){
		hist foreground(MC,"He3eta",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist background1(MC,"He3pi0pi0",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist background2(MC,"He3pi0pi0pi0",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		hist measured(DATA,"He3",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		value<double> norm(foreground.Total());
		if(norm.val()>1){
			foreground/=norm;			
			acceptance_fg[bin_num].varY()=norm/mc_norm_fg[bin_num].Y();
			background1/=value<double>(background1.Total());
			background2/=value<double>(background2.Total());
			
			Equation<DifferentialMutations<Parabolic>> fit([&measured,&foreground,&background1,&background2](const ParamSet&P)->double{
				return ChiSq(measured,(foreground*value<double>(P[0],0))
					+(background1*value<double>(P[1],0))
					+(background2*value<double>(P[2],0))
				,3);
			});
			fit.SetFilter([](const ParamSet&P)->bool{return (P[0]>=0)&&(P[1]>=0)&&(P[2]>=0);});
			fit.Init(200,make_shared<GenerateByGauss>()<<make_pair(norm.val()*2.0,norm.val()*2.0)<<make_pair(norm.val(),norm.val())<<make_pair(norm.val(),norm.val()),engine);
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
			PlotHist().Hist(measured,string("Data for bin where Q=")+to_string(mc_norm_fg[bin_num].X().val())+"+/-"+to_string(mc_norm_fg[bin_num].X().delta()))
				.Hist((foreground*value<double>(fit[0],0))+(background1*value<double>(fit[1],0))+(background2*value<double>(fit[2],0)),"Fit")
				.Hist(foreground*value<double>(fit[0],0),"^3He+eta").Hist(background1*value<double>(fit[1],0),"^3He+2pi^0").Hist(background2*value<double>(fit[2],0),"^3He+3pi^0")
			<<"set xlabel 'Missing mass, GeV'"<<"set ylabel 'Events count'";
			chi_sq[bin_num].varY()=value<double>(fit.Optimality(),0);
			events_fg[bin_num].varY()=value<double>(fit[0],errors[0]);
		}else{
			acceptance_fg[bin_num].varY()=value<double>(0,0);
			events_fg[bin_num].varY()=value<double>(0,0);
		}
	}
	PlotHist().Hist(chi_sq)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Fit chi^2, n.d.'";
	PlotHist().Hist(acceptance_fg,"^3He+eta")<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
	PlotHist().Hist(events_fg,"^3He+eta")<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	auto luminosity=events_fg/(acceptance_fg*sigmaHe3eta);
	PlotHist().Hist(luminosity)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integral luminosity, nb^{-1}'";
}
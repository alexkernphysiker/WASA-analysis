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
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
	vector<string> histpath_forward_reconstr={"Histograms","He3Forward_Reconstruction"};
	vector<string> reaction={"He3eta","He3pi0pi0pi0","He3pi0pi0"};
	vector<hist<double>> norm;
	for(const string& r:reaction)norm.push_back(Hist(MC,r,histpath_forward_reconstr,"0-Reference"));
	Plot<double>().Hist(norm[0],"Simulated events")
	.Hist(Hist(MC,reaction[0],histpath_forward_reconstr,"2-FPC"),"Forward tracks with signal in FPC")
	.Hist(Hist(MC,reaction[0],histpath_forward_reconstr,"3-AllCuts"),"identified as ^3He")
	<< "set key on"
	<< "set yrange [0:400000]"
	<< "set xlabel 'Q, MeV'"
	<< "set ylabel 'Events count'";

	vector<hist<double>> acceptance;
	for(const auto&h:norm)acceptance.push_back(hist<double>());
	hist<double> luminosity,bg_chi_sq,bg_ratio;
	RANDOM r_eng;
	for(size_t bin_num=0,bin_count=norm[0].size();bin_num<bin_count;bin_num++)
		if(norm[0][bin_num].X()>0.0){
			auto Q=norm[0][bin_num].X();
			string Qmsg="Q in ["+to_string(norm[0][bin_num].X().min())+":"+to_string(norm[0][bin_num].X().max())+"] MeV";
			auto transform=[](hist<double>&h){h=h.XRange(0.40,0.58);};

			hist<double> data=Hist(DATA,"",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
			transform(data);
			Plot<double>().Hist(data,"DATA "+Qmsg)
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'counts'"
			<< "set yrange [0:]";
		
			vector<hist<double>> theory;{
				for(size_t i=0;i<reaction.size();i++){
					hist<double> react_sim=Hist(MC,reaction[i],histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
					transform(react_sim);
					auto N=norm[i][bin_num].Y();
					acceptance[i] << point<value<double>>(Q,value<double>(react_sim.TotalSum().val())/N);
					react_sim/=N;
					theory.push_back(react_sim);
					react_sim/=2.0*react_sim[0].X().uncertainty();
					Plot<double>().Hist(theory[i],reaction[i]+" MC "+Qmsg)
					<< "set yrange [0:]"
					<< "set key on" 
					<< "set ylabel 'acceptance density, MeV^{-1}'";
				}
			}
			vector<LinearInterpolation<double>> reaction_funcs{theory[0].toLine(),theory[1].toLine(),theory[2].toLine()};
			Fit<DifferentialMutations<>,ChiSquare> fit(
				make_shared<FitPoints>(data),
				[&reaction_funcs](const ParamSet&X,const ParamSet&P){
					double res=0;
					for(size_t i=0;i<reaction_funcs.size();i++)
						res+=reaction_funcs[i](X[0])*P[i];
					return res;
				}
			);
			ParamSet pExit{	0.001,	0.001,	0.001	},
			pDelta{	0.01,	0.01,	0.01,	};
			fit.SetUncertaintyCalcDeltas(pDelta).SetFilter(make_shared<Above>()<<0.0<<0.0<<0.0);
			{
				auto count=data.TotalSum().val();
				fit.Init(300,make_shared<GenerateUniform>()<<make_pair(0.0,20.0*count)<<make_pair(0.0,20.0*count)<<make_pair(0.0,20.0*count),r_eng);
			}
			while((!fit.AbsoluteOptimalityExitCondition(0.000001))&&(!fit.ParametersDispersionExitCondition(pExit))){
				fit.Iterate(r_eng);
				cout<<fit.iteration_count()<<" iterations; "
				<<fit.Optimality()<<"<chi^2<"
				<<fit.Optimality(fit.PopulationSize()-1)
				<<"        \r";
			}
			const auto&P=fit.ParametersWithUncertainties();
			bg_ratio << point<value<double>>(Q,P[2]/P[1]);
			bg_chi_sq << point<value<double>>(Q,fit.Optimality()/(fit.Points()->size()-fit.ParamCount()));
			
			hist<double> FIT=theory[0]*P[0]+theory[1]*P[1]+theory[2]*P[2];
			Plot<double>()
			.Hist(data,"DATA "+Qmsg).Hist(FIT,"background")
			.Line(hist<double>(theory[0]*P[0]).toLine(),"^3He eta")
			.Line(hist<double>(theory[1]*P[1]).toLine(),"^3He3pi^0")
			.Line(hist<double>(theory[2]*P[2]).toLine(),"^3He2pi^0")
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'" 
			<< "set ylabel 'counts'"
			<< "set yrange [0:]";
			
			if(norm[0][bin_num].X()>5.0)
				luminosity << point<value<double>>(Q,(P[0]/func_value(he3eta_sigma().func(),Q))*double(trigger_he3_forward.scaling));
		}
	Plot<double>().Hist(bg_chi_sq) 
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'chi^2, n.d.'" 
	<< "set yrange [0:]";

	{//Plot acceptance
		Plot<double> plot;
		plot << "set key on" << "set yrange [0:1.0]";
		for(size_t i=0;i<reaction.size();i++)plot.Hist(acceptance[i],reaction[i]);
		plot << "set xlabel 'Q, MeV'" << "set ylabel 'Acceptance, n.d.'";
	}

	Plot<double>().Hist(bg_ratio/(acceptance[1]/acceptance[2])) 
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'sigma("+reaction[1]+")/sigma("+reaction[2]+"), n.d.'" 
	<< "set yrange [0:]";

	auto runs=PresentRuns("");
	Plot<double>().Hist(luminosity,to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs") 
	<< "set key on" << "set xlabel 'Q, MeV'" 
	<< "set ylabel 'Integral luminosity, nb^{-1}'" 
	<< "set yrange [0:100]";

	Plot<double>().Line(he3eta_sigma(),"Used in calculations")
	<< "set key on" << "set xlabel 'Q, MeV'" 
	<< "set ylabel 'sigma(^3He eta), nb'"<< "set yrange [0:600]";

	for(int i=2;i<=4;i++){
		auto phidistr=Hist(DATA,"",{"Histograms","He3Forward_Debug"},to_string(i)+"-PhiDistribution-AllBins").Scale(30);
		auto phidistr_mc=Hist(MC,"He3eta",{"Histograms","He3Forward_Debug"},to_string(i)+"-PhiDistribution-AllBins").Scale(30);
		Plot<double>().Hist(phidistr,"DATA")<<"set key on"
		<< "set xlabel 'Phi, deg'"<< "set ylabel 'events, count'"<< "set yrange [0:]";
		Plot<double>().Hist(phidistr_mc,"MC")<<"set key on"
		<< "set xlabel 'Phi, deg'"<< "set ylabel 'events, count'"<< "set yrange [0:]";
	}
}

// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Genetic/equation.h>
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
	vector<string> reaction={"He3eta","He3pi0pi0","He3pi0pi0pi0","He3pi0"};
	vector<hist<double>> norm;
	for(const string& r:reaction)norm.push_back(Hist(MC,r,histpath_forward_reconstr,"0-Reference"));

	vector<hist<double>> acceptance;
	for(const auto&h:norm)acceptance.push_back(hist<double>());
	hist<double> luminosity,bg_chi_sq,bg_ratio;
	RANDOM r_eng;
	for(size_t bin_num=0,bin_count=norm[0].size();bin_num<bin_count;bin_num++)
		if(norm[0][bin_num].X()>2.5){
			auto Q=norm[0][bin_num].X();
			string Qmsg="Q in ["+to_string(norm[0][bin_num].X().min())+":"+to_string(norm[0][bin_num].X().max())+"] MeV";
			auto transform=[](hist<double>&h){h=h.XRange(0.401,0.58);};

			hist<double> data=Hist(DATA,"",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
			transform(data);
			Plot<double> exp_plot;
			exp_plot.Hist(data,"DATA")
			<< "set key on"<< "set title '"+Qmsg+"'"
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'counts'"
			<< "set yrange [0:]";
		
			vector<hist<double>> theory;{
				Plot<double> th_plot;
				for(size_t i=0;i<reaction.size();i++){
					hist<double> react_sim=Hist(MC,reaction[i],histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
					transform(react_sim);
					auto N=norm[i][bin_num].Y();
					acceptance[i] << point<value<double>>(Q,value<double>(react_sim.TotalSum().val())/N);
					react_sim/=N;
					theory.push_back(react_sim);
					react_sim/=2.0*react_sim[0].X().uncertainty();
					th_plot.Line(react_sim.toLine(),reaction[i]);
				}
				th_plot<< "set yrange [0:]"<< "set key on"<< "set xlabel 'Missing mass, GeV'"
				<< "set title '"+Qmsg+"'"
				<< "set ylabel 'acceptance density, GeV^{-1}'";
				
			}
			vector<LinearInterpolation<double>> reaction_funcs{theory[0].toLine(),theory[1].toLine(),theory[2].toLine()};
			SearchMin<DifferentialMutations<ParabolicErrorEstimationFromChisq>> fit([&reaction_funcs,&data](const ParamSet&P){
				double res=0;
				for(size_t i=0,n=data.size();i<n;i++){
					value<double> exp_p=data[i].Y(),the_p=0;
					for(size_t j=0;j<reaction_funcs.size();j++)
						the_p+=reaction_funcs[j][i].Y()*P[j];
					res+=exp_p.NumCompare(the_p);
				}
				return res;
			});
			ParamSet pExit{	0.0001,	0.0001,	0.0001	},
			pDelta{	0.005,	0.005,	0.005	};
			fit.SetUncertaintyCalcDeltas(pDelta).SetFilter(make_shared<Above>()<<0.0<<0.0<<0.0);
			{
				auto count=data.TotalSum().val();
				fit.Init(300,make_shared<GenerateUniform>()<<make_pair(0.0,20.0*count)<<make_pair(0.0,20.0*count)<<make_pair(0.0,20.0*count),r_eng);
			}
			while((!fit.AbsoluteOptimalityExitCondition(0.0000001))&&(!fit.ParametersDispersionExitCondition(pExit))){
				fit.Iterate(r_eng);
				cout<<fit.iteration_count()<<" iterations; "
				<<fit.Optimality()<<"<chi^2<"
				<<fit.Optimality(fit.PopulationSize()-1)
				<<"        \r";
			}
			const auto&P=fit.ParametersWithUncertainties();
			bg_ratio << point<value<double>>(Q,
				(P[1]/acceptance[1].right().Y())
				/
				(P[2]/acceptance[2].right().Y())
			);
			bg_chi_sq << point<value<double>>(Q,fit.Optimality()/(data.size()-fit.ParamCount()));
			
			exp_plot.Line(hist<double>(theory[0]*P[0]+theory[1]*P[1]+theory[2]*P[2]).toLine(),"Total fit")
			.Line(hist<double>(theory[0]*P[0]).toLine(),"^3He eta")
			.Line(hist<double>(theory[1]*P[1]).toLine(),"^3He3pi^0")
			.Line(hist<double>(theory[2]*P[2]).toLine(),"^3He2pi^0");
			
			luminosity << point<value<double>>(Q,(
				(P[0]/acceptance[0].right().Y())
				/func_value(he3eta_sigma().func(),Q)
			)*double(trigger_he3_forward.scaling));
		}
	Plot<double>().Hist(bg_chi_sq) 
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'chi^2/d, n.d.'" 
	<< "set yrange [0:]";

	{//Plot acceptance
		Plot<double> plot;
		plot << "set key on" 
		<< "set title 'Acceptance'"
		<< "set yrange [0:1.0]"
		<< "set xlabel 'Q, MeV'" 
		<< "set ylabel 'Acceptance, n.d.'";
		for(size_t i=0;i<reaction.size();i++)plot.Hist(acceptance[i],reaction[i]);
	}

	Plot<double>().Hist(bg_ratio) 
	<< "set title 'Background reactions'"
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'sigma("+reaction[1]+")/sigma("+reaction[2]+"), n.d.'" 
	<< "set yrange [0:]";

	auto runs=PresentRuns("");
	Plot<double>().Hist(luminosity,to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs") 
	<< "set title 'Integral luminosity'"
	<< "set key on" << "set xlabel 'Q, MeV'" 
	<< "set ylabel 'Integral luminosity, nb^{-1}'" 
	<< "set yrange [0:]";

	Plot<double>().Line(he3eta_sigma(),"Used in calculations")
	<< "set title '"+reaction[0]+"'"
	<< "set key on" << "set xlabel 'Q, MeV'" 
	<< "set ylabel 'sigma(^3He eta), nb'"<< "set yrange [0:600]";
}

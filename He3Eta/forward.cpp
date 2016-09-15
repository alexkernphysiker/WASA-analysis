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
	vector<string> reaction={"He3eta","He3pi0pi0","He3pi0pi0pi0"};
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
		if(norm[0][bin_num].X()>17.5){
			auto Q=norm[0][bin_num].X();
			string Qmsg="Q in ["+to_string(norm[0][bin_num].X().min())+":"+to_string(norm[0][bin_num].X().max())+"] MeV";
			auto transform=[](hist<double>&h){h=h.XRange(0.35,0.75);};

			hist<double> data=Hist(DATA,"",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
			transform(data);
			Plot<double>()
			.Hist(data,"DATA "+Qmsg)
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'counts'"
			<< "set xrange [0.4:0.6]"
			<< "set yrange [-200:2500]";
		
			vector<hist<double>> theory;{
				for(size_t i=0;i<reaction.size();i++){
					hist<double> react_sim=Hist(MC,reaction[i],histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
					transform(react_sim);
					auto N=norm[i][bin_num].Y();
					acceptance[i] << point<value<double>>(Q,value<double>(react_sim.TotalSum().val())/N);
					theory.push_back(react_sim/N);
					Plot<double>().Hist(theory[i],reaction[i]+" MC "+Qmsg)
					<< "set xrange [0.4:0.6]"
					<< "set yrange [0:]"
					<< "set key on" 
					<< "set ylabel 'acceptance density, channel^{-1}'";
				}
			}
			vector<LinearInterpolation<double>> bg_funcs{theory[1].toLine(),theory[2].toLine()};
			Fit<DifferentialMutations<>,ChiSquareWithXError> bg_fit(
				make_shared<FitPoints>(data.XRange(0.450,0.530)),
				[&bg_funcs](const ParamSet&X,const ParamSet&P){
					double res=0;
					for(size_t i=0;i<bg_funcs.size();i++)res+=bg_funcs[i](X[0])*P[i];
					return res;
				}
			);
			bg_fit.SetUncertaintyCalcDeltas({0.01,0.01}).SetFilter(make_shared<Above>()<<0.0<<0.0);
			{
				auto count=data.TotalSum().val();
				bg_fit.Init(100,make_shared<GenerateUniform>()
					<<make_pair(0.0,20.0*count)
					<<make_pair(0.0,20.0*count)
					,r_eng
				);
			}
			while(!bg_fit.AbsoluteOptimalityExitCondition(0.0000001))
				bg_fit.Iterate(r_eng);
			bg_ratio << point<value<double>>(Q,
				bg_fit.ParametersWithUncertainties()[0]
				/
				bg_fit.ParametersWithUncertainties()[1]
			);
			
			SortedPoints<double> BG_displ(theory[1].toLine()*bg_fit.Parameters()[0]+theory[2].toLine()*bg_fit.Parameters()[1]);
			Plot<double>()
			.Hist(bg_fit.Points()->Hist1(0),"cut DATA "+Qmsg)
			.Line(BG_displ,"fit")
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'counts'"
			<< "set xrange [0.4:0.6]"
			<< "set yrange [-200:2500]";
			
			bg_chi_sq << point<value<double>>(Q,bg_fit.Optimality()/(bg_fit.Points()->size()-bg_fit.ParamCount()));
			
			hist<double> BG=
				theory[1]*bg_fit.ParametersWithUncertainties()[0]+
				theory[2]*bg_fit.ParametersWithUncertainties()[1];
			Plot<double>()
			.Hist(data,"DATA "+Qmsg).Hist(BG,"background")
			.Line(hist<double>(theory[1]*bg_fit.ParametersWithUncertainties()[0]).toLine(),"^3He2pi^0")
			.Line(hist<double>(theory[2]*bg_fit.ParametersWithUncertainties()[1]).toLine(),"^3He3pi^0")
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'" 
			<< "set ylabel 'counts'"
			<< "set xrange [0.4:0.6]"
			<< "set yrange [-200:2500]";
			
			hist<double> FG=(data-BG).XRange(0.40,0.70);
			Plot<double>().Object("0*x title ''")
			.Hist(FG,"DATA-background "+Qmsg)
			<< "set key on" 
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'counts'"
			<< "set xrange [0.4:0.6]"
			<< "set yrange [-200:2500]";
			
			FG=FG.XRange(0.525,0.560);
			value<double> L=FG.TotalSum()/theory[0].TotalSum().val();
			Plot<double>().Object("0*x title ''")
			.Hist(FG,"DATA-background "+Qmsg)
			.Line(hist<double>(theory[0]*L).toLine(),"^3He+eta MC*N")
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'counts'"
			<< "set xrange [0.4:0.6]"
			<< "set yrange [-200:2500]";
			luminosity << point<value<double>>(Q,
				L*double(trigger_he3_forward.scaling)
				/
				func_value(he3eta_sigma().func(),Q)
			);
		}
	Plot<double>().Hist(bg_chi_sq) 
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'chi^2, n.d.'" 
	<< "set yrange [0:10]";
			
	{//Plot acceptance
		Plot<double> plot;
		plot << "set key on" << "set yrange [0:0.7]";
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
	<< "set yrange [0:80]";
	
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

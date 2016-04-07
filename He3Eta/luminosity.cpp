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
#include <Kinematics/particles.h>
#include <Kinematics/reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
const Reaction&main_reaction(){
	static Reaction main_react(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	return main_react;
}
int main(){
	LinearInterpolation<double> sigmaHe3eta{
		//http://arxiv.org/pdf/nucl-ex/0701072v1
		point<double>(-0.5,0.0),
		point<double>(0.0,100.0),
		point<double>(0.5,380.0),
		point<double>(1.5,400.0),
		point<double>(6.0,390.0),
		point<double>(12.0,380.0),
		//Extrapolating
		point<double>(24.0,360.0),
		point<double>(36.0,340.0)
	};
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
	auto Q2P=LinearInterpolation<double>(SortedPoints<double>([](double p){return main_reaction().P2Q(p);},ChainWithStep(0.0,0.001,3.0)).Transponate());
	auto Q2E=LinearInterpolation<double>(SortedPoints<double>([](double e){return main_reaction().E2Q(e);},ChainWithStep(0.0,0.001,3.0)).Transponate());
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
	for(const auto&h:norm)acceptance.push_back(h.CloneEmptyBins());
	SortedPoints<value<double>> luminosity,bg_chi_sq;
	RANDOM r_eng;
	for(size_t bin_num=3,bin_count=norm[0].size();bin_num<bin_count;bin_num++){
		auto transform=[](hist<double>&h){h=h.Scale(3).XRange(0.45,0.6);};
		Plotter::Instance() << "unset yrange" << "unset xrange";
		hist<double> data=Hist(DATA,"He3",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
		transform(data);
		Plot<double>().Hist(data,"DATA")
		<< "set key on"
		<< "set xlabel 'Missing mass, GeV'"
		<< "set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
		<< "set xrange [0.45:0.6]"
		<< "set yrange [0:]";
		
		vector<hist<double>> theory;
		for(size_t i=0;i<reaction.size();i++){
			hist<double> react_sim=Hist(MC,reaction[i],histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
			transform(react_sim);
			auto N=value<double>(react_sim.Total());
			acceptance[i].Bin(bin_num).varY()=N/norm[i][bin_num].Y();
			theory.push_back(react_sim/norm[i][bin_num].Y());
		}
		Plot<double>().Hist(theory[0],"^3He+eta (MC)") << "set xrange [0.45:0.6]"<< "set key on";
		if(norm[0][bin_num].X()>20.0){
			vector<LinearInterpolation<double>> bg_funcs{theory[1].Line(),theory[2].Line()};
			Fit<DifferentialMutations<>,ChiSquareWithXError> bg_fit(
				make_shared<FitPoints>(data.XRange(0.47,0.56).XExclude(0.532,0.555)),
				[&bg_funcs](const ParamSet&X,const ParamSet&P){
					double res=0;
					for(size_t i=0;i<bg_funcs.size();i++)res+=bg_funcs[i](X[0])*P[i];
					return res;
				}
			);
			bg_fit.SetUncertaintyCalcDeltas({0.01,0.01}).SetFilter(make_shared<Above>()<<0.0<<0.0);
			auto count=data.Total()*20.0;
			bg_fit.Init(100,make_shared<GenerateUniform>()<<make_pair(0.0,count)<<make_pair(0.0,count),r_eng);
			while(!bg_fit.AbsoluteOptimalityExitCondition(0.0000001))
				bg_fit.Iterate(r_eng);
			
			SortedPoints<double> BG_displ(theory[1].Line()*bg_fit[0]+theory[2].Line()*bg_fit[1]);
			Plot<double>().Hist(bg_fit.Points()->Hist1(0),"cut DATA").Line(BG_displ,"fit")
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<< "set xrange [0.45:0.6]"
			<< "set yrange [0:]";
			
			bg_chi_sq << point<value<double>>(norm[0][bin_num].X(),bg_fit.Optimality());
			
			hist<double> BG=
				theory[1]*bg_fit.ParametersWithUncertainties()[0]+
				theory[2]*bg_fit.ParametersWithUncertainties()[1];
			Plot<double>().Hist(data,"all DATA").Hist(BG,"background MC")
			.Line(hist<double>(theory[1]*bg_fit.ParametersWithUncertainties()[0]).Line(),"^3He2pi^0")
			.Line(hist<double>(theory[2]*bg_fit.ParametersWithUncertainties()[1]).Line(),"^3He3pi^0")
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'" 
			<< "set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<< "set xrange [0.45:0.6]"
			<< "set yrange [0:]";
			
			hist<double> FG=(data-BG).XRange(0.49,0.60);
			Plot<double>().Hist(FG,"DATA-background").Object("0*x title ''")
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'"<<"set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<< "set xrange [0.45:0.6]"
			<< "set yrange [-300:2000]";
			
			FG=FG.XRange(0.530,0.555);
			value<double> L=FG.TotalSum()/theory[0].Total();
			Plot<double>().Hist(FG,"DATA-background").Hist(theory[0]*L,"^3He+eta MC*N")
			<< "set key on"
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<< "set xrange [0.45:0.6]"
			<< "set yrange [-300:2000]";
			luminosity << point<value<double>>(
				norm[0][bin_num].X(),
				L*value<double>(trigger_he3_forward.scaling)
				/
				func_value(sigmaHe3eta.func(),norm[0][bin_num].X())
			);
		}
	}
	Plot<double>().Hist(bg_chi_sq) 
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'BG fit chi^2, n.d.'" 
	<< "set yrange [0:]";
	{//Plot acceptance
		Plot<double> plot;
		plot << "set key on" << "set yrange [0:0.7]";
		for(size_t i=0;i<reaction.size();i++)plot.Hist(acceptance[i],reaction[i]);
		plot << "set xlabel 'Q, MeV'" << "set ylabel 'Acceptance, n.d.'";
	}
	auto runs=PresentRuns("He3");
	Plot<double>().Hist(luminosity,to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs") 
	<< "set key on" << "set xlabel 'Q, MeV'" 
	<< "set ylabel 'Integral luminosity, nb^{-1}'" 
	<< "set yrange [0:100]";
}

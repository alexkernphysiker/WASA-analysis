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
		<<"set yrange [0:400000]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	Plotter::Instance()<<"unset yrange";
	vector<hist<double>> acceptance;
	for(const auto&h:norm)acceptance.push_back(h.CloneEmptyBins());
	SortedPoints<value<double>> luminosity,bg_chi_sq;
	RANDOM r_eng;
	for(size_t bin_num=3,bin_count=norm[0].size();bin_num<bin_count;bin_num++){
		auto transform=[](hist<double>&h){h=h.Scale(2).XRange(0.4,0.7);};
		Plotter::Instance()<<"unset yrange"<<"unset xrange";
		hist<double> measured=Hist(DATA,"He3",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
		transform(measured);
		Plot<double>().Hist(measured,"DATA")
		<< "set xlabel 'Missing mass, GeV'"
		<< "set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
		<< "set xrange [0.4:0.6]"
		<< "set yrange [0:]";
		
		vector<hist<double>> theory;
		for(size_t i=0;i<reaction.size();i++){
			hist<double> react_sim=Hist(MC,reaction[i],histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
			transform(react_sim);
			auto N=value<double>(react_sim.Total());
			acceptance[i].Bin(bin_num).varY()=N/norm[i][bin_num].Y();
			theory.push_back(react_sim/norm[i][bin_num].Y());
		}
		Plot<double>().Hist(theory[0]);
		auto eta_exists=theory[0].YRange(0.005,+INFINITY);
		if((eta_exists.size()>0)&&(norm[0][bin_num].X()>17.0)){
			vector<LinearInterpolation<double>> bg_funcs{theory[1].Line(),theory[2].Line()};
			Fit<DifferentialMutations<>,ChiSquareWithXError> bg_fit(
				make_shared<FitPoints>(
					measured.XRange(0.50,0.60)
					.XExclude(eta_exists.left().X(),eta_exists.right().X())
				),
				[&bg_funcs](const ParamSet&X,const ParamSet&P){
					double res=0;
					for(size_t i=0;i<bg_funcs.size();i++)res+=bg_funcs[i](X[0])*P[i];
					return res;
				}
			);
			bg_fit.SetUncertaintyCalcDeltas({0.01,0.01}).SetFilter(make_shared<Above>()<<0.0<<0.0);
			auto count=measured.Total()*20.0;
			bg_fit.Init(100,make_shared<GenerateUniform>()<<make_pair(0.0,count)<<make_pair(0.0,count),r_eng);
			while(!bg_fit.AbsoluteOptimalityExitCondition(0.0000001))
				bg_fit.Iterate(r_eng);
			
			SortedPoints<double> BG_displ(theory[1].Line()*bg_fit[0]+theory[2].Line()*bg_fit[1]);
			Plot<double>().Hist(bg_fit.Points()->Hist1(0),"cut DATA").Line(BG_displ,"fit")
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<< "set xrange [0.4:0.6]"
			<< "set yrange [0:]";
			
			bg_chi_sq << point<value<double>>(norm[0][bin_num].X(),bg_fit.Optimality());
			
			hist<double> BG=theory[1]*bg_fit.ParametersWithUncertainties()[0]+theory[2]*bg_fit.ParametersWithUncertainties()[1];
			Plot<double>().Hist(measured,"all DATA").Hist(BG,"background")
			.Line(hist<double>(theory[1]*bg_fit.ParametersWithUncertainties()[0]).Line(),"^3He2pi^0")
			<< "set xlabel 'Missing mass, GeV'" 
			<< "set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"<<"set yrange [0:]";
			
			Plotter::Instance() << "unset yrange" << "unset xrange";
			hist<double> FG=(measured-BG).XRange(0.50,0.60);
			Plot<double>().Hist(FG,"substracted").Object("0*x title ''")
			<< "set xlabel 'Missing mass, GeV'"<<"set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<< "set xrange [0.5:0.6]"
			<< "set yrange [-300:1500]";
			
			FG=FG.XRange(0.53,0.56).YRange(0.0,+INFINITY);
			value<double> L=FG.TotalSum()/theory[0].Total();
			Plot<double>().Hist(FG,"substracted").Hist(theory[0].XRange(0.50,0.60)*L,"MC")
			<< "set xlabel 'Missing mass, GeV'"
			<< "set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<< "set xrange [0.5:0.6]"
			<< "set yrange [-300:2000]";
			luminosity << point<value<double>>(
				norm[0][bin_num].X(),
				L*value<double>(trigger_he3_forward.scaling)
				/
				func_value(sigmaHe3eta.func(),norm[0][bin_num].X())
			);
		}
	}
	Plot<double>().Hist(bg_chi_sq) << "set xlabel 'Q, MeV'" << "set ylabel 'BG fit chi^2, n.d.'" << "set yrange [0:]";
	{//Plot acceptance
		Plot<double> plot;
		plot<<"set yrange [0:0.7]";
		for(size_t i=0;i<reaction.size();i++)plot.Hist(acceptance[i],reaction[i]);
		plot << "set xlabel 'Q, MeV'" << "set ylabel 'Acceptance, n.d.'";
	}
	Plot<double>().Hist(luminosity) << "set xlabel 'Q, MeV'" << "set ylabel 'Integral luminosity, a.u.'" << "set yrange [0:]";
}

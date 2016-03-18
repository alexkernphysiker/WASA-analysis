// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/functions.h>
#include <math_h/error.h>
#include <math_h/interpolate.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Kinematics/particles.h>
#include <Kinematics/reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace GnuplotWrap;
const Reaction&main_reaction(){
	static Reaction main_react(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	return main_react;
}
Plot<double> cs_plot;
value<double> sigmaHe3eta(const value<double>&Q){
	static LinearInterpolation<double> sigma;
	if(sigma.size()==0){
		sigma
		//http://arxiv.org/pdf/nucl-ex/0701072v1
		<<point<double>(-0.5,0.0)
		<<point<double>(0.0,100.0)
		<<point<double>(0.5,380.0)
		<<point<double>(1.5,400.0)
		<<point<double>(6.0,390.0)
		<<point<double>(12.0,380.0)
		//Extrapolating
		<<point<double>(24.0,360.0)
		<<point<double>(36.0,340.0);
		cs_plot.Line(sigma,"3He eta");
		cs_plot<<"set xlabel 'Q, MeV'"<<"set ylabel 'cross section, mb'";
	}
	return sigma.func()*Q;
}
value<double> sigmaHe3pi0pi0(const value<double>&Q){
	static LinearInterpolation<double> sigma;
	if(sigma.size()==0){
		sigma
		<<point<double>(-0.5,2800.0)
		<<point<double>(32.0,2800.0);
		cs_plot.Line(sigma,"3He 2pi0");
	}
	return sigma.func()*Q;
}
value<double> sigmaHe3pi0pi0pi0(const value<double>&Q){
	static LinearInterpolation<double> sigma;
	if(sigma.size()==0){
		sigma
		//10.1140/epja/i2010-10981-3
		<<point<double>(-0.5,115.0)
		<<point<double>(32.0,115.0);
		cs_plot.Line(sigma,"3He 3pi0");
	}
	return sigma.func()*Q;
}
int main(){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
	auto Q2P=LinearInterpolation<double>(SortedPoints<double>([](double p){return main_reaction().P2Q(p);},ChainWithStep(0.0,0.001,3.0)).Transponate());
	auto Q2E=LinearInterpolation<double>(SortedPoints<double>([](double e){return main_reaction().E2Q(e);},ChainWithStep(0.0,0.001,3.0)).Transponate());
	vector<string> histpath_forward={"Histograms","He3Forward_Reconstruction"};
	vector<string> reaction={"He3eta","He3pi0pi0","He3pi0pi0pi0"};
	vector<hist<double>::Func> cross_section={sigmaHe3eta,sigmaHe3pi0pi0,sigmaHe3pi0pi0pi0};
	vector<hist<double>> norm;
	for(const string& r:reaction)norm.push_back(Hist(MC,r,histpath_forward,"0-Reference"));
	Plot<double>().Hist(norm[0],"All events")
		.Hist(Hist(MC,reaction[0],histpath_forward,"1-AllTracks"),"All tracks in forward")
		.Hist(Hist(MC,reaction[0],histpath_forward,"2-FPC"),"Signal in FPC")
		.Hist(Hist(MC,reaction[0],histpath_forward,"3-AllCuts"),"^3He")
		<<"set yrange [0:400000]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	Plotter::Instance()<<"unset yrange";
	vector<hist<double>> acceptance;
	for(const auto&h:norm)acceptance.push_back(h.CloneEmptyBins());
	SortedPoints<value<double>> luminosity;
	for(size_t bin_num=5,bin_count=norm[0].size();bin_num<bin_count;bin_num++){
		Plot<double> mc_plot;
		hist<double> theory;
		Plotter::Instance()<<"unset yrange"<<"unset xrange";
		for(size_t i=0;i<reaction.size();i++){
			hist<double> react_sim=Hist(MC,reaction[i],histpath_forward,string("MissingMass-Bin-")+to_string(bin_num)).XRange(0.4,0.6);
			auto N=value<double>(react_sim.Total());
			acceptance[i].Bin(bin_num).varY()=N/norm[i][bin_num].Y();
			react_sim/=norm[i][bin_num].Y();
			react_sim*=cross_section[i];
			if(theory.size()==0)
				theory=react_sim;
			else
				theory+=react_sim;
			mc_plot.Hist(react_sim,reaction[i]);
		}
		mc_plot.Line(theory.Line(),"Sum")<<"set xlabel 'Missing mass, GeV'"<<"set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"<<"set yrange [0:]";

		hist<double> measured=Hist(DATA,"He3",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num)).XRange(0.4,0.6);
		auto K=value<double>(measured.Total())/theory.TotalSum();
		Plot<double>().Hist(measured,"DATA").Hist(theory*K,"Simulation")
			<<"set xlabel 'Missing mass, GeV'"<<"set ylabel 'counts (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<<"set yrange [0:]";
		luminosity<<point<value<double>>(norm[0][bin_num].X(),K*value<double>(trigger_he3_forward.scaling,0));
	}
	{//Plot acceptance
		Plot<double> plot;
		plot<<"set yrange [0:0.7]";
		for(size_t i=0;i<reaction.size();i++)plot.Hist(acceptance[i],reaction[i]);
		plot<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
	}
	Plotter::Instance()<<"unset yrange";
	Plot<double>().Hist(luminosity)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integral luminosity, a.u.'"<<"set yrange [0:]";
}

// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/functions.h>
#include <math_h/error.h>
#include <Genetic/fit.h>
#include <Genetic/equation.h>
#include <Genetic/initialconditions.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Kinematics/particles.h>
#include <Kinematics/reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace GnuplotWrap;
using namespace Genetic;
RANDOM random_engine;
const Reaction&main_reaction(){
	static Reaction main_react(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	return main_react;
}
Plot<double> cs_plot;
double sigmaHe3eta(const double Q){
	static LinearInterpolation<double> sigma;
	if(sigma.size()==0){
		sigma
		//http://arxiv.org/pdf/nucl-ex/0701072v1
		<<make_pair(-0.5,0.0)
		<<make_pair(0.0,100.0)
		<<make_pair(0.5,380.0)
		<<make_pair(1.5,400.0)
		<<make_pair(6.0,390.0)
		<<make_pair(12.0,380.0)
		//Extrapolating
		<<make_pair(24.0,360.0)
		<<make_pair(36.0,340.0);
		cs_plot.Line(sigma,"3He eta");
		cs_plot<<"set xlabel 'Q, MeV'"<<"set ylabel 'cross section, mb'";
	}
	return sigma(Q);
}
double sigmaHe3pi0pi0(const double Q){
	static LinearInterpolation<double> sigma;
	if(sigma.size()==0){
		sigma
		//Proposal
		<<make_pair(-0.5,2800.0)
		<<make_pair(32.0,2800.0);
		cs_plot.Line(sigma,"3He 2pi0");
	}
	return sigma(Q);
}
double sigmaHe3pi0pi0pi0(const double Q){
	static LinearInterpolation<double> sigma;
	if(sigma.size()==0){
		sigma
		//10.1140/epja/i2010-10981-3
		<<make_pair(-0.5,115.0)
		<<make_pair(32.0,115.0);
		cs_plot.Line(sigma,"3He 3pi0");
	}
	return sigma(Q);
}
int main(){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
	auto Q2P=LinearInterpolation<double>([](double p){return main_reaction().P2Q(p);},ChainWithStep(0.0,0.001,3.0)).Transponate();
	auto Q2E=LinearInterpolation<double>([](double e){return main_reaction().E2Q(e);},ChainWithStep(0.0,0.001,3.0)).Transponate();
	vector<string> histpath_forward={"Histograms","He3Forward_Reconstruction"};
	vector<string> reaction={"He3eta","He3pi0pi0","He3pi0pi0pi0"};
	vector<function<double(double)>> cross_section={sigmaHe3eta,sigmaHe3pi0pi0,sigmaHe3pi0pi0pi0};
	vector<hist<double>> norm;
	for(const string& r:reaction)norm.push_back(Hist(MC,r,histpath_forward,"0-Reference"));
	Plot<double>().Hist(norm[0],"All events")
		.Hist(Hist(MC,reaction[0],histpath_forward,"1-AllTracks"),"Forward Tracks")
		.Hist(Hist(MC,reaction[0],histpath_forward,"2-FPC"),"Signal in FPC")
		.Hist(Hist(MC,reaction[0],histpath_forward,"3-AllCuts"),"E_{dep} cuts")
		.Hist(Hist(MC,reaction[0],histpath_forward,"5-Kinematic cut"),"Kinematic cuts")
		<<"set yrange [0:400000]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	Plotter::Instance()<<"unset yrange";
	vector<hist<double>> acceptance;
	for(const auto&h:norm)acceptance.push_back(h.CloneEmptyBins());
	vector<point<double>> luminosity;
	for(size_t bin_num=5,bin_count=norm[0].size();bin_num<bin_count;bin_num++){
		vector<hist<double>> simulations;
		Plotter::Instance()<<"unset yrange"<<"unset xrange";
		{
			Plot<double> mc_plot;
			for(size_t i=0;i<reaction.size();i++){
				auto react_sim=Hist(MC,reaction[i],histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
				auto N=value<double>(react_sim.Total());
				acceptance[i].Bin(bin_num).varY()=N/norm[i][bin_num].Y();
				react_sim/=norm[i][bin_num].Y();
				simulations.push_back(react_sim);
				mc_plot.Hist(react_sim,reaction[i]);
			}
			mc_plot<<"set xlabel 'Missing mass, GeV'"
				<<"set ylabel 'a.u (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
				<<"set yrange [0:]"<<"set logscale y";
		}
		Plotter::Instance()<<"unset yrange"<<"unset xrange"<<"unset logscale";
		auto measured=Hist(DATA,"He3",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		SearchMin<DifferentialMutations<Parabolic>> find_parameters(
			[&simulations,&measured](const ParamSet&P)->double{
				auto sum=measured.CloneEmptyBins();
				for(size_t i=0,n=P.size();i<n;i++)
					sum+=simulations[i]*value<double>(P[i],0);
				return measured.ChiSq_only_y_error(sum,P.size());
			}
		);
		find_parameters.SetFilter([](const ParamSet&P)->bool{
			for(double p:P)if(p<0.0)return false;
			return true;
		});
		find_parameters.Init(60,make_shared<GenerateByGauss>()<<make_pair(1000.0,1000.0)<<make_pair(1000.0,1000.0)<<make_pair(1000.0,1000.0),random_engine);
		while(!find_parameters.AbsoluteOptimalityExitCondition(0.0000001))
			find_parameters.Iterate(random_engine);
		
		Plot<double>().Hist(measured,"DATA")
		.Hist(
			simulations[0]*value<double>(find_parameters[0],0)+
			simulations[1]*value<double>(find_parameters[1],0)+
			simulations[2]*value<double>(find_parameters[2],0)
			,"fit"
		)
		<<"set xlabel 'Missing mass, GeV'"<<"set ylabel 'counts (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
		<<"set yrange [0:]";
		luminosity.push_back(point<double>(
			norm[0][bin_num].X(),
			value<double>(find_parameters[0],find_parameters.GetParamParabolicError(0.1,0))/
			value<double>(cross_section[0](norm[0][bin_num].X().val()),0)
		));
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
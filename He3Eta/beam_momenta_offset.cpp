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
const Reaction&main_reaction(){
	static Reaction main_react(Particle::p(),Particle::d(),{Particle::he3(),Particle::eta()});
	return main_react;
}
int main(){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"beam_momenta");
	auto Q2P=LinearInterpolation<double>([](double p){return main_reaction().P2Q(p);},ChainWithStep(0.0,0.001,3.0)).Transponate();
	auto Q2E=LinearInterpolation<double>([](double e){return main_reaction().E2Q(e);},ChainWithStep(0.0,0.001,3.0)).Transponate();
	vector<string> histpath_forward={"Histograms","He3Forward_Reconstruction"};
	string he3eta="He3eta";
	Plot<double> theory_plot;
	theory_plot<<"set xlabel 'E_{kin}, GeV'"<<"set ylabel 'theta, deg'";
	auto QBins=Hist(MC,he3eta,histpath_forward,"0-Reference");
	for(const auto&P:QBins)
		theory_plot.Line(
			LinearInterpolation<double>(
				[&P,&Q2P](double E)->double{
					return main_reaction().PbEr2Theta(Q2P(P.X().val()/1000.0),E)*180./3.1415926;
				},
			       ChainWithStep(0.2,0.001,0.4)
			),
		to_string(P.X().val())
		);
	RANDOM engine;
	vector<point<double>> offs_vertex,offs_mc,offs_data;
	for(size_t bin_num=9,bin_count=QBins.size();bin_num<bin_count;bin_num++){
		auto x=QBins[bin_num].X();
		auto do_fit=[&engine,&x,&Q2P](const hist2d<double>&kin_)->point<double>{
			auto points=make_shared<FitPoints>();
			double max=0;
			kin_.FullCycle([&max](point3d<double>&&P){
				if(max<P.Z().val())max=P.Z().val();
			});
			kin_.FullCycle([max,&points](point3d<double>&&P){
				if((P.Z().val()>(2.0*max/3.0))&&(P.X().val()>0.25)&&(P.X().val()<0.35))
					points<<Point({P.X().val()},P.Y().val(),P.Z().val());
			});
			Fit<DifferentialMutations<>,SumWeightedSquareDiff> fit(
				points,
				[&x,&Q2P](const ParamSet&E,const ParamSet&P){
					return main_reaction().PbEr2Theta(Q2P((x.val()+P[0])/1000.0),E[0])*180./3.1415926;
				}
			);
			fit.Init(50,make_shared<GenerateUniform>()<<make_pair(-10.0,10.0),engine);
			while(!fit.AbsoluteOptimalityExitCondition(0.000001))
				fit.Iterate(engine);
			Plot<double> kin_plot;
			kin_plot.Points(points->Hist1(0).Line()).Line(
				LinearInterpolation<double>(
					[&fit](double E)->double{return fit({E});},
					ChainWithStep(0.2,0.001,0.4)
				)
			);
			return point<double>(x,value<double>(fit[0],0));
		};
		auto kin_v=Hist2d(MC,he3eta,histpath_forward,string("Kinematic-vertex-Bin-")+to_string(bin_num)).Scale(4,4);
		PlotHist2d<double>(sp2).Distr(kin_v)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'";
		offs_vertex.push_back(do_fit(kin_v));
		auto kin_mc=Hist2d(MC,he3eta,histpath_forward,string("Kinematic-before-cut-Bin-")+to_string(bin_num)).Scale(4,4);
		PlotHist2d<double>(sp2).Distr(kin_mc)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'";
		offs_mc.push_back(do_fit(kin_mc));
		auto kin_data=Hist2d(DATA,"He3",histpath_forward,string("Kinematic-before-cut-Bin-")+to_string(bin_num)).Scale(4,4);
		PlotHist2d<double>(sp2).Distr(kin_data);
		offs_data.push_back(do_fit(kin_data));
	}
	Plot<double>().Hist(offs_vertex,"Vertex").Hist(offs_mc,"WMC").Hist(offs_data,"Data");
	return 0;
}
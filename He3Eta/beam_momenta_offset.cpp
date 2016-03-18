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
	auto Q2P=LinearInterpolation<double>(SortedPoints<double>([](double p){return main_reaction().P2Q(p);},ChainWithStep(0.0,0.005,3.0)).Transponate());
	auto Q2E=LinearInterpolation<double>(SortedPoints<double>([](double e){return main_reaction().E2Q(e);},ChainWithStep(0.0,0.005,3.0)).Transponate());
	vector<string> histpath_forward={"Histograms","He3Forward_Reconstruction"};
	string he3eta="He3eta";
	RANDOM engine;
	SortedPoints<value<double>> offs_mc,offs_data;
	auto QBins=Hist(MC,he3eta,histpath_forward,"0-Reference");
	for(size_t bin_num=10,bin_count=QBins.size();bin_num<bin_count;bin_num++){
		auto Q=QBins[bin_num].X();
		double p=Q2P(Q.val()/1000.0);
		auto do_fit=[&engine,&Q,&p](const hist2d<double>&kin_)->point<value<double>>{
			auto points=make_shared<FitPoints>();
			double max=0;
			kin_.FullCycle([&max](const point3d<value<double>>&P){
				if(max<P.Z().val())max=P.Z().val();
			});
			kin_.FullCycle([max,&points](const point3d<value<double>>&P){
				if((P.Z().val()>(2.0*max/3.0))&&(P.X().val()>0.28)&&(P.X().val()<0.35))
					points<<Point({P.X().val()},P.Y().val(),P.Z().val());
			});
			Fit<DifferentialMutations<>,SumWeightedSquareDiff> fit(
				points,
				[p](const ParamSet&E,const ParamSet&P){
					return main_reaction().PbEr2Theta(p+P[0],E[0])*180./PI();
				}
			);
			fit.Init(30,make_shared<GenerateUniform>()<<make_pair(-0.01,0.01),engine);
			while(!fit.AbsoluteOptimalityExitCondition(0.00001))
				fit.Iterate(engine);
			Plot<double> kin_plot;
			SortedPoints<double> pts;
			for(const auto&pt:points->Hist1(0))pts<<point<double>(pt.X().val(),pt.Y().val());
			kin_plot.Points(pts).Line(
				SortedPoints<double>(
					[&fit](double E)->double{return fit({E});},
					ChainWithStep(0.2,0.005,0.4)
				)
			);
			return point<value<double>>(Q,value<double>(fit[0],0));
		};
		auto kin_v=Hist2d(MC,he3eta,histpath_forward,string("Kinematic-vertex-Bin-")+to_string(bin_num)).Scale(4,4);
		PlotHist2d<double>(sp2).Distr(kin_v)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'";
		auto kin_mc=Hist2d(MC,he3eta,histpath_forward,string("Kinematic-reconstructed-Bin-")+to_string(bin_num)).Scale(4,4);
		PlotHist2d<double>(sp2).Distr(kin_mc)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'";
		offs_mc<<do_fit(kin_mc);
		auto kin_data=Hist2d(DATA,"He3",histpath_forward,string("Kinematic-reconstructed-Bin-")+to_string(bin_num)).Scale(4,4);
		PlotHist2d<double>(sp2).Distr(kin_data);
		offs_data<<do_fit(kin_data);
	}
	Plot<double>().Hist(offs_mc,"WMC").Hist(offs_data,"Data")<<"set yrange [0:.006]"
		<<"set xlabel 'Q, MeV'"<<"set ylabel 'delta, GeV'";
	return 0;
}
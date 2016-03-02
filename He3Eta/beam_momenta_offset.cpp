// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/functions.h>
#include <math_h/error.h>
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
	for(size_t bin_num=9,bin_count=QBins.size();bin_num<bin_count;bin_num++){
		auto x=QBins[bin_num].X();
		auto do_plot=[&x,&Q2P](const hist2d<double>&kin_){
			vector<pair<double,double>> points;
			double max=0;kin_.FullCycle([&max](point3d<double>&&P){
				if(max<P.Z().val())
					max=P.Z().val();
			});
			kin_.FullCycle([max,&points](point3d<double>&&P){
				if(P.Z().val()>(2.0*max/3.0))
					points.push_back(make_pair(P.X().val(),P.Y().val()));
			});
			Plot<double> kin_plot;
			kin_plot.Points(points);
			for(double offset=0.000;offset<0.010;offset+=0.001)
				kin_plot.Line(
					LinearInterpolation<double>(
						[&offset,&x,&Q2P](double E)->double{
							return main_reaction().PbEr2Theta(Q2P(x.val()/1000.0+offset),E)*180./3.1415926;
						},
				 ChainWithStep(0.2,0.001,0.4)
					),
		  to_string(offset)
				);
		};
		auto kin_v=Hist2d(MC,he3eta,histpath_forward,string("Kinematic-vertex-Bin-")+to_string(bin_num));
		kin_v=kin_v.Scale(4,4);
		PlotHist2d<double>(sp2).Distr(kin_v)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'";
		do_plot(kin_v);
		auto kin_mc=Hist2d(MC,he3eta,histpath_forward,string("Kinematic-before-cut-Bin-")+to_string(bin_num));
		kin_mc=kin_mc.Scale(4,4);
		PlotHist2d<double>(sp2).Distr(kin_mc)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'";
		do_plot(kin_mc);
		auto kin_data=Hist2d(DATA,"He3",histpath_forward,string("Kinematic-before-cut-Bin-")+to_string(bin_num));
		kin_data=kin_data.Scale(4,4);
		PlotHist2d<double>(sp2).Distr(kin_data);
		do_plot(kin_data);
	}
	return 0;
}
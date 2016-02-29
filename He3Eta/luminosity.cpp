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
		cs_plot.Line(sigma,"3He+eta");
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
		cs_plot.Line(sigma,"3He+2pi0");
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
		cs_plot.Line(sigma,"3He+3pi0");
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
	{
		Plot<double> kin_plot;
		for(const auto&P:norm[0])
			kin_plot.Line(
				LinearInterpolation<double>(
					[&P,&Q2P](double E)->double{
						return main_reaction().PbEr2Theta(Q2P(P.X().val()/1000.0),E)*180./3.1415926;
					},
					ChainWithStep(0.2,0.001,0.4)
				),
				to_string(P.X().val())
			);
	}
	for(size_t bin_num=5,bin_count=norm[0].size();bin_num<bin_count;bin_num++){
		Plot<double> mc_plot1;
		Plot<double> mc_plot;
		hist<double> theory;
		Plotter::Instance()<<"unset yrange"<<"unset xrange";
		{
			auto x=norm[0][bin_num].X();
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
			auto kin_v=Hist2d(MC,reaction[0],histpath_forward,string("Kinematic-vertex-Bin-")+to_string(bin_num));
			kin_v=kin_v.Scale(4,4);
			PlotHist2d<double>(sp2).Distr(kin_v)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'";
			do_plot(kin_v);
			auto kin_mc=Hist2d(MC,reaction[0],histpath_forward,string("Kinematic-before-cut-Bin-")+to_string(bin_num));
			kin_mc=kin_mc.Scale(4,4);
			PlotHist2d<double>(sp2).Distr(kin_mc)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'";
			do_plot(kin_mc);
			auto kin_data=Hist2d(DATA,"He3",histpath_forward,string("Kinematic-before-cut-Bin-")+to_string(bin_num));
			kin_data=kin_data.Scale(4,4);
			PlotHist2d<double>(sp2).Distr(kin_data);
			do_plot(kin_data);
		}
		for(size_t i=0;i<reaction.size();i++){
			auto react_sim=Hist(MC,reaction[i],histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
			mc_plot1.Hist(react_sim,reaction[i]);
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

		auto measured=Hist(DATA,"He3",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
		auto K=value<double>(measured.Total())/theory.TotalSum();
		Plot<double>().Hist(measured,"DATA").Line((theory*K).Line(),"Simulation")
			<<"set xlabel 'Missing mass, GeV'"<<"set ylabel 'counts (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<<"set yrange [0:]";
		luminosity.push_back(point<double>(norm[0][bin_num].X(),K*value<double>(trigger_he3_forward.scaling,0)));
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

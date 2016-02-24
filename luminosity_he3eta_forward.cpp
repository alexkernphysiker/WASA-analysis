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
#include <str_get.h>
#include <gethist.h>
#include <particles.h>
#include <reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace Genetic;
using namespace GnuplotWrap;
int main(int,char**){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
	vector<string> histpath_forward={"Histograms","He3Forward_Reconstruction"};
	vector<string> reaction={"He3eta","He3pi0pi0","He3pi0pi0pi0"};
	vector<function<double(double)>> cross_section={
		[](double q)->double{return sigmaHe3eta(PBeam_He3eta(q));},
		[](double q)->double{return sigmaHe3pi0pi0(PBeam_He3eta(q));},
		[](double q)->double{return sigmaHe3pi0pi0pi0(PBeam_He3eta(q));}
	};
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
	vector<point<double>> luminosity,mc_offs,data_offs;
	{
		Plot<double> kin_plot;
		for(const auto&P:norm[0])
			kin_plot.Line(
				LinearInterpolation<double>(
					[&P](double E)->double{
						return Ekin2Theta_He3eta(E,PBeam_He3eta(P.X().val()/1000.0));
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
		if(bin_num>8){
			auto x=norm[0][bin_num].X();
			RANDOM r;
			auto do_fit=[&x,&r](const hist2d<double>&kin_)->point<double>{
				auto points=make_shared<FitPoints>();
				double max=0;kin_.FullCycle([&max](point3d<double>&&P){
					if(max<P.Z().val())
						max=P.Z().val();
				});
				kin_.FullCycle([max,&points](point3d<double>&&P){
					if(P.Z().val()>(max/2.0))
						points<<Point({P.X().val()},P.Y().val(),P.Z().val());
				});
				Fit<DifferentialMutations<>,SumWeightedSquareDiff> fit(
					points,
					[&x](const ParamSet&E,const ParamSet&Offs)->double{
						return Ekin2Theta_He3eta(E[0],PBeam_He3eta(x.val()/1000.0)-Offs[0]);
					}
				);
				fit.Init(50,make_shared<GenerateByGauss>()<<make_pair(0,0.01),r);
				while(!fit.RelativeOptimalityExitCondition(0.00001))fit.Iterate(r);
				Plot<double>().Points(points->Hist1(0).Line())
				.Line(LinearInterpolation<double>(
					[&fit](double e)->double{return fit({e});},
					ChainWithStep(0.25,0.001,0.35)
				));
				return point<double>(x,value<double>(fit[0],0));
			};
			auto kin_mc=Hist2d(MC,reaction[0],histpath_forward,string("Kinematic-before-cut-Bin-")+to_string(bin_num));
			kin_mc=kin_mc.Scale(4,4);
			PlotHist2d<double>(sp2).Distr(kin_mc)<<"set xlabel 'E_k, GeV'"<<"set ylabel 'theta, deg'";
			mc_offs.push_back(do_fit(kin_mc));
			auto kin_data=Hist2d(DATA,"He3",histpath_forward,string("Kinematic-before-cut-Bin-")+to_string(bin_num));
			kin_data=kin_data.Scale(4,4);
			PlotHist2d<double>(sp2).Distr(kin_data);
			data_offs.push_back(do_fit(kin_data));
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
		luminosity.push_back(point<double>(norm[0][bin_num].X(),K));
	}
	Plot<double>().Hist(mc_offs,"MC").Hist(data_offs,"DATA");
	{//Plot acceptance
		Plot<double> plot;
		plot<<"set yrange [0:0.7]";
		for(size_t i=0;i<reaction.size();i++)plot.Hist(acceptance[i],reaction[i]);
		plot<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
	}
	Plotter::Instance()<<"unset yrange";
	Plot<double>().Hist(luminosity)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integral luminosity, a.u.'"<<"set yrange [0:]";
}
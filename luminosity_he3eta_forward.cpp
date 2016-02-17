// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/functions.h>
#include <math_h/error.h>
#include <str_get.h>
#include <gethist.h>
#include <particles.h>
#include <reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace GnuplotWrap;
int main(int,char**){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
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
	for(size_t bin_num=6,bin_count=norm[0].size()-1;bin_num<bin_count;bin_num++){
		Plot<double> mc_plot1;
		Plot<double> mc_plot;
		hist<double> theory;
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
	{//Plot acceptance
		Plot<double> plot;
		plot<<"set yrange [0:0.7]";
		for(size_t i=0;i<reaction.size();i++)plot.Hist(acceptance[i],reaction[i]);
		plot<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
	}
	Plot<double>().Hist(hist<double>(luminosity))<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integral luminosity, a.u.'"<<"unset yrange"<<"set yrange [0:]";
}
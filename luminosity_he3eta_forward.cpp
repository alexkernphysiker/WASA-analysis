// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <math_h/functions.h>
#include <math_h/error.h>
#include <Genetic/fit.h>
#include <Genetic/equation.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
#include <str_get.h>
#include <gethist.h>
#include <particles.h>
#include <reactions.h>
using namespace std;
using namespace ROOT_data;
using namespace GnuplotWrap;
using namespace Genetic;
int main(int,char**){
	RANDOM engine;
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
	vector<string> histpath_forward={"Histograms","He3Forward_Reconstruction"};
	vector<string> reaction={"He3eta","He3pi0pi0"};
	vector<function<double(double)>> cross_section={sigmaHe3eta,sigmaHe3pi0pi0};
	vector<hist> norm;
	for(const string& r:reaction)norm.push_back(hist(MC,r,histpath_forward,"0-Reference"));
	vector<hist> acceptance;
	for(const auto&h:norm)acceptance.push_back(h.CloneEmptyBins());
	auto luminosity=norm[0].CloneEmptyBins();
	for(size_t bin_num=6;bin_num<norm[0].count();bin_num++){
		vector<hist> simulated;
		for(const string r:reaction)simulated.push_back(hist(MC,r,histpath_forward,string("MissingMass-Bin-")+to_string(bin_num)));
		vector<value<double>> mc_event_count;
		PlotHist mc_plot;
		for(size_t i=0;i<reaction.size();i++){
			mc_event_count.push_back(value<double>(simulated[i].Total()));
			acceptance[i][bin_num].varY()=mc_event_count[i]/norm[i][bin_num].Y();
			simulated[i]/=norm[i][bin_num].Y();
			simulated[i]*=cross_section[i];
			mc_plot.Hist(simulated[i],reaction[i]);
		}
		mc_plot<<"set xlabel 'Missing mass, GeV'"<<"set ylabel 'counts (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"<<"set yrange [0:]";
		if(mc_event_count[0].val()>1){
			hist measured(DATA,"He3",histpath_forward,string("MissingMass-Bin-")+to_string(bin_num));
			hist theory=measured.CloneEmptyBins();
			for(const hist&H:simulated)theory+=H;
			luminosity[bin_num].varY()=value<double>(measured.Total())/value<double>(theory.Total());
			PlotHist().Hist(measured,"DATA").Hist(theory*luminosity[bin_num].Y(),"Simulation")
			<<"set xlabel 'Missing mass, GeV'"<<"set ylabel 'counts (Q="+to_string(norm[0][bin_num].X().val())+" MeV)'"
			<<"set yrange [0:]";
		}
	}
	PlotHist().Hist(luminosity)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integral luminosity, a.u.'"<<"set yrange [0:]";
	{//Plot acceptance
		PlotHist plot;
		plot<<"set yrange [0:0.7]";
		for(size_t i=0;i<reaction.size();i++)plot.Hist(acceptance[i],reaction[i]);
		plot<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'";
	}
	{// debug messaging
		hist n1(MC,reaction[0],histpath_forward,"1-AllTracks");
		hist n2(MC,reaction[0],histpath_forward,"2-FPC");
		hist n3(MC,reaction[0],histpath_forward,"3-AllCuts");
		hist n4(MC,reaction[0],histpath_forward,"4-Reconstructed");
		hist n5(MC,reaction[0],histpath_forward,"5-Kinematic cut");
		PlotHist().Hist(norm[0],"All events").Hist(n1,"Forward Tracks")
		.Hist(n2,"Signal in FPC").Hist(n3,"E_{dep} cuts").Hist(n5,"Kinematic cuts")
		<<"set yrange [0:400000]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	}
}
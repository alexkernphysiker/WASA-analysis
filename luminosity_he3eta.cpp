// this file is distributed under 
// GPL v 3.0 license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <fit.h>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
#include <equation.h>
#include <str_get.h>
#include <gethist.h>
#include <theory.h>
#include <histfit.h>
using namespace std;
using namespace Genetic;
RANDOM engine;
#define Q_range 15.0,30.0
#define MissingMass_range 0.5,0.6
int main(int,char**){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta");
	hist mc_norm(MC,"He3eta",{"Histograms","Reconstruction"},"Reference");
	hist acceptance(MC,"He3eta",{"Histograms","Reconstruction"},"Additional");
	{
		hist mc_filtered1(MC,"He3eta",{"Histograms","Reconstruction"},"Theta_reconstruction_correct");
		hist mc_filtered2(MC,"He3eta",{"Histograms","Reconstruction"},"Reconstructed");
		PlotHist().Hist("All MC events",mc_norm).Hist("FPC",mc_filtered1)
			.Hist("Reconstructed",mc_filtered2).Hist("Preselected",acceptance)
			<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	}
	(acceptance/=mc_norm).Cut(Q_range);
	PlotHist().Hist("Acceptance",acceptance)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'"<<"set nokey";
	hist luminosity=acceptance.CloneEmptyBins();
	for(auto&qBin:luminosity){
		int index=int(qBin.x*1000.0);
		vector<hist> MChists={
			hist(MC,"He3eta",{"Histograms","MissingMass"},to_string(index)),
			hist(MC,"He3pi0pi0",{"Histograms","MissingMass"},to_string(index)),
			hist(MC,"He3pi0pi0pi0",{"Histograms","MissingMass"},to_string(index))
		};
		hist DATAhist(DATA,"He3",{"Histograms","MissingMass"},to_string(index));
		for(hist&H:MChists)H.Cut(MissingMass_range);
		DATAhist.Cut(MissingMass_range);
		string suffix=string("Q=")+to_string(qBin.x)+"MeV";
		PlotHist().HistWLine(string("MCHe3eta")+suffix,MChists[0])
			.HistWLine(string("MCHe3 2pi0")+suffix,MChists[1])
			.HistWLine(string("MCHe3 3pi0")+suffix,MChists[2])
			<<"set xlabel 'MM, GeV'"<<"set ylabel 'Counts'";
		hist result=FitHistByHists(DATAhist,MChists,engine,{&(qBin.y)},{&(qBin.dy)});
		PlotHist().Hist(string("DataHe3")+suffix,DATAhist)
			.HistWLine(string("Fit")+suffix,result)
			<<"set xlabel 'MM, GeV'"<<"set ylabel 'Counts'";
	}
	const double MC_events_count=5000000;
	PlotHist().Hist("He3eta true events in data",luminosity*=MC_events_count)
		<<"set xlabel 'Q, MeV'"<<"set ylabel 'True events, counts'"<<"set nokey";
	PlotHist lumplot;
	lumplot.Hist("analysed",luminosity/=sigmaHe3eta);
	lumplot.Hist("estimated",luminosity/=PresentRunsAmountRatio("He3"));
	lumplot<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integrated luminocity, 1/nb'";
}
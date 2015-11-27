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
int main(int,char**){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta");
	hist mc_norm(false,"He3eta",{"Histograms","Reconstruction"},"Reference");
	hist acceptance(false,"He3eta",{"Histograms","Reconstruction"},"Additional");
	{
		hist mc_filtered1(false,"He3eta",{"Histograms","Reconstruction"},"Theta_reconstruction_correct");
		hist mc_filtered2(false,"He3eta",{"Histograms","Reconstruction"},"Reconstructed");
		PlotHist().Hist("All MC events",mc_norm).Hist("FPC",mc_filtered1)
		.Hist("Reconstructed",mc_filtered2).Hist("Preselected",acceptance);
	}
	(acceptance/=mc_norm).Cut(15,30);
	PlotHist().Hist("Acceptance",acceptance);
	hist luminosity=acceptance.CloneEmptyBins();
	
	for(auto&qBin:luminosity){
		int index=int(qBin.x*1000);
		vector<hist> MC={
			hist(false,"He3eta",{"Histograms","MissingMass"},to_string(index)),
			hist(false,"He3pi0pi0",{"Histograms","MissingMass"},to_string(index)),
			hist(false,"He3pi0pi0pi0",{"Histograms","MissingMass"},to_string(index))
		};
		hist data(true,"He3",{"Histograms","MissingMass"},to_string(index));
		for(hist&H:MC)
			H.Cut(0.5,0.6);
		data.Cut(0.5,0.6);
		hist result=FitHistByHists(
			data,MC,engine,
			{&(qBin.y)},{&(qBin.dy)}
		);
		string suffix=string("Q=")+to_string(qBin.x)+"MeV";
		PlotHist().HistWLine(string("MCHe3eta")+suffix,MC[0])
			.HistWLine(string("MCHe3 2pi0")+suffix,MC[1])
			.HistWLine(string("MCHe3 3pi0")+suffix,MC[2]);
		PlotHist().Hist(string("DataHe3")+suffix,data)
			.HistWLine(string("Fit")+suffix,result);
	}
	PlotHist().Hist("He3eta true events in data",luminosity);
	PlotHist lumplot;
	lumplot.Hist("Integral luminosity (analysed)",luminosity/=sigmaHe3eta);
	lumplot.Hist("Integral luminosity (estimated)",luminosity/=PresentRunsAmountRatio("He3"));
}
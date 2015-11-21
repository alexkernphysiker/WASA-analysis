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
	
	hist mc_norm(false,"He3eta",{"Histograms","Cuts"},"Reference"),
	mc_filtered1(false,"He3eta",{"Histograms","Cuts"},"IsInFPC"),
	mc_filtered2(false,"He3eta",{"Histograms","Cuts"},"Edep_cuts"),
	acceptance(false,"He3eta",{"Histograms","Cuts"},"Reconstructed");
	PlotHist().Hist("All MC events",mc_norm).Hist("FPC",mc_filtered1)
	.Hist("E_{dep} cuts",mc_filtered2).Hist("Reconstructed",acceptance);
	//ToDo: refactor luminosity
	(acceptance/=mc_norm).Cut(1.61,1.65);
	PlotHist().Hist("Acceptance",acceptance);
	hist luminosity=acceptance.CloneEmptyBins();
	
	for(auto&BeamMomentaBin:luminosity){
		int index=int(BeamMomentaBin.x*1000);
		vector<hist> MC={
			hist(false,"He3eta",{"Histograms","MissingMass"},to_string(index)),
			hist(false,"He3pi0",{"Histograms","MissingMass"},to_string(index)),
			hist(false,"He3pi0pi0",{"Histograms","MissingMass"},to_string(index)),
			hist(false,"He3pi0pi0pi0",{"Histograms","MissingMass"},to_string(index))
		};
		hist data(true,"He3",{"Histograms","MissingMass"},to_string(index));
		data.Cut(0.5,0.6);for(hist&H:MC)H.Cut(0.5,0.6);
		hist result=FitHistByHists(BeamMomentaBin,data,MC,engine);
		string suffix=string(" P=")+to_string(BeamMomentaBin.x)+"GeVc";
		{hist bg1=MC[1],bg2=MC[2];
			PlotHist().HistWLine(string("MCHe3eta")+suffix,MC[0])
				.Hist(string("MCHe3 2pi0")+suffix,bg1*=5)
				.Hist(string("MCHe3 3pi0")+suffix,bg2*=5);
		}
		PlotHist().Hist(string("DataHe3")+suffix,data)
			.HistWLine(string("Fit")+suffix,result);
	}
	PlotHist().Hist("He3eta true events in data",luminosity*=1000000.0);
	PlotHist lumplot;
	lumplot.Hist("Integral luminosity (analysed)",luminosity/=sigmaHe3eta);
	lumplot.Hist("Integral luminosity (estimated)",luminosity/=PresentRunsAmountRatio("He3"));
}
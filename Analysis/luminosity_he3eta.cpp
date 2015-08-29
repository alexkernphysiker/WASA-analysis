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
#include "gethist.h"
#include "fit_hist2hist.h"
using namespace std;
using namespace Genetic;
RANDOM engine;
double sigma(double p_beam){
	return 400;//nb
}
void AnalyseMMSpectra(hist::point&BeamMomentaBin,hist&data,vector<hist>&MC){
	for(auto p:data)if(p.x>=0.542)BeamMomentaBin.y+=p.y;
	BeamMomentaBin.dy=sqrt(BeamMomentaBin.y);
	if(BeamMomentaBin.dy<1)BeamMomentaBin.dy=1;
}
int main(int,char**){
#include "env.cc"
	Plotter::Instance().SetOutput(outpath,"he3eta");
	
	hist mc_norm(false,"He3eta",{"Histograms","Cuts"},"Reference"),
	mc_filtered1(false,"He3eta",{"Histograms","Cuts"},"EdepCuts"),
	mc_filtered2(false,"He3eta",{"Histograms","Cuts"},"ThetaCut"),
	acceptance(false,"He3eta",{"Histograms","Cuts"},"Reconstructed");
	PlotHist().Hist("All MC events",mc_norm).Hist("EdepCut",mc_filtered1)
		.Hist("ThetaCut",mc_filtered2).Hist("Reconstructed",acceptance);
	PlotHist().Hist("Acceptance",acceptance/=mc_norm);

	hist luminocity=acceptance.CloneEmptyBins();
	for(auto&BeamMomentaBin:luminocity){
		string suffix=to_string(int(BeamMomentaBin.x*1000));
		vector<hist> MC={
			hist(false,"He3eta",{"Histograms","MissingMass"},static_cast<string&&>(suffix)),
			hist(false,"He3pi0pi0",{"Histograms","MissingMass"},static_cast<string&&>(suffix)),
			hist(false,"He3pi0pi0pi0",{"Histograms","MissingMass"},static_cast<string&&>(suffix))
		};
		PlotHist()
			.Hist(string("MCHe3eta")+suffix,static_cast<hist&&>(MC[0]))
			.Hist(string("MCHe3 2pi0")+suffix,static_cast<hist&&>(MC[1]))
			.Hist(string("MCHe3 3pi0")+suffix,static_cast<hist&&>(MC[2]));
		hist data(true,"He3eta",{"Histograms","MissingMass"},static_cast<string&&>(suffix));
		PlotHist().Hist(string("DataHe3eta")+suffix,data);
		AnalyseMMSpectra(BeamMomentaBin,data,MC);
	}
	PlotHist().Hist("He3eta events in data",luminocity)
		.Hist("He3eta true events",luminocity/=acceptance);
	PlotHist().Hist("Integral luminocity(analysed)",luminocity/=sigma)
		.Hist("Integral luminocity(estimated)",luminocity/=PresentRunsAmountRatio("He3eta"));
}
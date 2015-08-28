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
	vector<string> kin_path={"Histograms","MissingMass"};
	vector<string> norm_path={"Histograms","FRH1"};
	
	hist mc_norm(false,"He3eta",static_right(norm_path),"Reference"),
		mc_filtered1(false,"He3eta",static_right(norm_path),"SP2cut"),
		mc_filtered2(false,"He3eta",static_right(norm_path),"ThetaCut"),
		acceptance(false,"He3eta",static_right(norm_path),"Reconstructed");
	PlotHist().Hist("All MC events",static_right(mc_norm))
		.Hist("EdepCut",static_right(mc_filtered1))
		.Hist("ThetaCut",static_right(mc_filtered2))
		.Hist("Reconstructed",static_right(acceptance));
	acceptance/=mc_norm;
	PlotHist().Hist("Acceptance",static_right(acceptance));
	hist luminocity=acceptance.CloneEmptyBins();
	for(auto&BeamMomentaBin:luminocity){
		string suffix=to_string(int(BeamMomentaBin.x*1000));
		vector<hist> MC={
			hist(false,"He3eta",static_right(kin_path),suffix),
			hist(false,"He3pi0pi0",static_right(kin_path),suffix),
			hist(false,"He3pi0pi0pi0",static_right(kin_path),suffix)
		};
		PlotHist()
			.Hist(string("MCHe3eta")+suffix,static_cast<hist&&>(MC[0]))
			.Hist(string("MCHe3 2pi0")+suffix,static_cast<hist&&>(MC[1]))
			.Hist(string("MCHe3 3pi0")+suffix,static_cast<hist&&>(MC[2]));
		hist data(true,"He3eta",static_right(kin_path),suffix);
		PlotHist().Hist(string("DataHe3eta")+suffix,static_right(data));
		AnalyseMMSpectra(BeamMomentaBin,data,MC);
	}
	{
		PlotHist eventcnt;
		eventcnt.Hist("He3eta events in data",static_right(luminocity));
		luminocity/=acceptance;
		eventcnt.Hist("He3eta true events",static_right(luminocity));
	}{
		PlotHist lum_plot;
		luminocity/=sigma;
		lum_plot.Hist("Integral luminocity(analysed)",static_right(luminocity));
		luminocity/=PresentRunsAmountRatio("He3eta");
		lum_plot.Hist("Integral luminocity(estimated)",static_right(luminocity));
	}
}
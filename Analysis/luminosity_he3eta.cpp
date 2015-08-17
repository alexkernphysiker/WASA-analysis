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
	vector<string> kin_path={"Histograms","Kinematics"};
	vector<string> norm_path={"Histograms","EventsCount"};
	
	hist mc_norm(false,"He3eta",static_right(norm_path),"AllEventsOnPBeam"),
		mc_filtered(false,"He3eta",static_right(norm_path),"FilteredEventsOnPBeam"),
		acceptance(false,"He3eta",static_right(norm_path),"ReconstructedEventsOnPBeam");
	PlotHist().Hist("All MC events",static_right(mc_norm))
		.Hist("Preselected",static_right(mc_filtered))
		.Hist("Reconstructed",static_right(acceptance));
	acceptance/=mc_norm;
	PlotHist().Hist("Acceptance",static_right(acceptance));
	string missingmass="MissingMass";
	hist luminocity=acceptance.CloneEmptyBins();
	for(hist::point&BeamMomentaBin:luminocity)if(BeamMomentaBin.y>=0.05){
		string suffix=to_string(int(BeamMomentaBin.x*1000));
		vector<hist> MC={
			hist(false,"He3eta",static_right(kin_path),missingmass+suffix),
			hist(false,"He3pi0pi0",static_right(kin_path),missingmass+suffix),
			hist(false,"He3pi0pi0pi0",static_right(kin_path),missingmass+suffix)
		};
		PlotHist()
			.Hist(string("MCHe3eta")+suffix,static_cast<hist&&>(MC[0]))
			.Hist(string("MCHe3 2pi0")+suffix,static_cast<hist&&>(MC[1]))
			.Hist(string("MCHe3 3pi0")+suffix,static_cast<hist&&>(MC[2]));
		hist data(true,"He3eta",static_right(kin_path),missingmass+suffix);
		PlotHist().Hist(string("DataHe3eta")+suffix,static_right(data));
		AnalyseMMSpectra(BeamMomentaBin,data,MC);
	}
	{
		PlotHist eventcnt;
		eventcnt.Hist("He3eta events in data",static_right(luminocity));
		luminocity/=acceptance;
		eventcnt.Hist("He3eta true events",static_right(luminocity));
	}
	luminocity/=sigma;
	PlotHist().Hist("Luminocity",static_right(luminocity));
}
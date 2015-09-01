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
using namespace std;
using namespace Genetic;
RANDOM engine;
double sigmaHe3eta(double p_beam){
	return 400;//nb
}
void AnalyseMMSpectra(hist::point&BeamMomentaBin,const hist&data,const vector<hist>&MC){
	string suffix=string(" P=")+to_string(BeamMomentaBin.x)+"GeVc";
	{hist bg1=MC[1],bg2=MC[2];
		PlotHist().HistWLine(string("MCHe3eta")+suffix,MC[0])
			.Hist(string("MCHe3 2pi0")+suffix,bg1*=5)
			.Hist(string("MCHe3 3pi0")+suffix,bg2*=5);
	}
	PlotHist fitplot;
	fitplot.Hist(string("DataHe3")+suffix,data);
	for(auto p:data)if(p.x>=0.542)BeamMomentaBin.y+=p.y;
	BeamMomentaBin.dy=sqrt(BeamMomentaBin.y);
	if(BeamMomentaBin.dy<1)BeamMomentaBin.dy=1;
}
int main(int,char**){
#include "env.cc"
	Plotter::Instance().SetOutput(outpath,"he3eta");
	
	hist mc_norm(false,"He3eta",{"Histograms","Cuts"},"Reference"),
	mc_filtered1(false,"He3eta",{"Histograms","Cuts"},"Edep_cuts"),
	mc_filtered2(false,"He3eta",{"Histograms","Cuts"},"ThetaCut"),
	acceptance(false,"He3eta",{"Histograms","Cuts"},"Reconstructed");
	PlotHist().Hist("All MC events",mc_norm).Hist("EdepCut",mc_filtered1)
		.Hist("ThetaCut",mc_filtered2).Hist("Reconstructed",acceptance);
	
	PlotHist().Hist("Acceptance",acceptance/=mc_norm);
	hist luminocity=acceptance.CloneEmptyBins();
	
	for(auto&BeamMomentaBin:luminocity){
		int index=int(BeamMomentaBin.x*1000);
		vector<hist> MC={
			hist(false,"He3eta",{"Histograms","MissingMass"},to_string(index)),
			hist(false,"He3pi0pi0",{"Histograms","MissingMass"},to_string(index)),
			hist(false,"He3pi0pi0pi0",{"Histograms","MissingMass"},to_string(index))
		};
		//Read data
		hist data(true,"He3",{"Histograms","MissingMass"},to_string(index));
		AnalyseMMSpectra(BeamMomentaBin,data,MC);
	}
	
	Plotter::Instance()<<"set xrange [1.58:*]";
	PlotHist eventsplot;
	eventsplot.Hist("He3eta events in data",luminocity);
	eventsplot.Hist("He3eta true events",luminocity/=acceptance);
	PlotHist lumplot;
	lumplot.Hist("Integral luminocity(analysed)",luminocity/=sigmaHe3eta);
	lumplot.Hist("Integral luminocity(estimated)",luminocity/=PresentRunsAmountRatio("He3"));
}
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
	auto histsum=[&data,&MC](const ParamSet&P){
		hist sum=data.CloneEmptyBins();
		for(size_t i=0;i<MC.size();i++){
			hist tmp=MC[i];
			tmp*=P[i];
			sum+=tmp;
		}
		return sum;
	};
	SearchMin<DifferentialMutations<Parabolic>> fit([&data,&MC,histsum](const ParamSet&P){
		hist bg=histsum(P);
		return data.HowClose(bg);
	});
	fit.SetFilter([](const ParamSet&P){
		bool res=true;
		for(double p:P)res&=(p>0);
		return res;
	});
	auto init=make_shared<GenerateByGauss>();
	for(auto H:MC)init<<make_pair(1,1);
	fit.Init(MC.size()*20,init,engine);
	while(!fit.AbsoluteOptimalityExitCondition(0.0000001))
		fit.Iterate(engine);
	hist fithist=histsum(fit.Parameters());
	fitplot.HistWLine(string("Fit")+suffix,fithist);
	BeamMomentaBin.y=fit[0];
	BeamMomentaBin.dy=fit.GetParamParabolicError(0.0000001,0);
}
int main(int,char**){
#include "env.cc"
	Plotter::Instance().SetOutput(outpath,"he3eta");
	
	hist mc_norm(false,"He3eta",{"Histograms","Cuts"},"Reference"),
	mc_filtered1(false,"He3eta",{"Histograms","Cuts"},"IsInFPC"),
	mc_filtered2(false,"He3eta",{"Histograms","Cuts"},"Edep_cuts"),
	acceptance(false,"He3eta",{"Histograms","Cuts"},"Reconstructed");
	PlotHist().Hist("All MC events",mc_norm).Hist("FPC",mc_filtered1)
		.Hist("E_{dep} cuts",mc_filtered2).Hist("Reconstructed",acceptance);
	
	(acceptance/=mc_norm).Cut(1.61,1.65);
	PlotHist().Hist("Acceptance",acceptance);
	hist luminosity=acceptance.CloneEmptyBins();
	
	for(auto&BeamMomentaBin:luminosity){
		int index=int(BeamMomentaBin.x*1000);
		vector<hist> MC={
			hist(false,"He3eta",{"Histograms","MissingMass"},to_string(index)),
			hist(false,"He3pi0pi0",{"Histograms","MissingMass"},to_string(index)),
			hist(false,"He3pi0pi0pi0",{"Histograms","MissingMass"},to_string(index))
		};
		//Read data
		hist data(true,"He3",{"Histograms","MissingMass"},to_string(index));
		data.Cut(0.5,0.6);for(hist&H:MC)H.Cut(0.5,0.6);
		AnalyseMMSpectra(BeamMomentaBin,data,MC);
	}
	PlotHist().Hist("He3eta true events in data",luminosity*=1000000.0);
	PlotHist lumplot;
	lumplot.Hist("Integral luminosity (analysed)",luminosity/=sigmaHe3eta);
	lumplot.Hist("Integral luminosity (estimated)",luminosity/=PresentRunsAmountRatio("He3"));
}
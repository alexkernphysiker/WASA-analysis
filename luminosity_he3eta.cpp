// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <functions.h>
#include <fit.h>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
#include <equation.h>
#include <str_get.h>
#include <gethist.h>
#include <theory.h>
#include "phys_constants.h"
using namespace std;
using namespace Genetic;
RANDOM engine;
#define Q_range 15.0,30.0
#define MissingMass_range 0.4,0.6
const double MC_events_count=5000000;
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
		hist MC_He3eta(MC,"He3eta",{"Histograms","MissingMass"},to_string(index)),
			MC_He3pi0pi0(MC,"He3pi0pi0",{"Histograms","MissingMass"},to_string(index)),
			MC_He3pi0pi0pi0(MC,"He3pi0pi0pi0",{"Histograms","MissingMass"},to_string(index)),
			DATAhist(DATA,"He3",{"Histograms","MissingMass"},to_string(index));
		MC_He3eta.Cut(MissingMass_range);
		DATAhist.Cut(MissingMass_range);
		string suffix=string("Q=")+to_string(qBin.x)+"MeV";
		PlotHist()
			.HistWLine(string("MCHe3eta")+suffix,MC_He3eta)
			.HistWLine(string("MCHe3 2pi0")+suffix,MC_He3pi0pi0)
			.HistWLine(string("MCHe3 3pi0")+suffix,MC_He3pi0pi0pi0)
			<<"set xlabel 'MM, GeV'"<<"set ylabel 'Counts'";
		printf("%f MeV \n",qBin.x);
		const size_t pop_size=50;
		
		typedef Mul<Par<0>,Func3<Gaussian,Arg<0>,Par<1>,Par<2>>> Peak;
		
		FitFunction<DifferentialMutations<>,Peak,ChiSquareWithXError> fitMC(make_shared<FitPoints>()<<MC_He3eta);
#define init make_pair(m_eta,0.01)<<make_pair(0.01,0.01)
		fitMC.SetFilter(make_shared<Above>()<<0<<0<<0)
			.Init(pop_size,make_shared<GenerateByGauss>()<<make_pair(100,100)<<init,engine);
		while(!fitMC.AbsoluteOptimalityExitCondition(0.0001))
			fitMC.Iterate(engine);
		PlotFit1D<decltype(fitMC)>()
			.Fit("mc_fit_"+suffix,"mc_"+suffix,fitMC,0.001)
			<<"set xlabel 'MM, GeV'"<<"set ylabel 'counts'";

		Fit<DifferentialMutations<>,ChiSquareWithXError> fitdata(
			make_shared<FitPoints>()<<DATAhist,
			[&MC_He3pi0pi0,&MC_He3pi0pi0pi0](const ParamSet&X,const ParamSet&P){
				static Peak peak;
				return peak(X,P)+P[3]*MC_He3pi0pi0(X[0]).y+P[4]*MC_He3pi0pi0pi0(X[0]).y;
			}
		);
		fitdata.SetFilter(make_shared<Above>()<<0<<0<<0<<0<<0)
			.Init(pop_size,make_shared<GenerateByGauss>()<<make_pair(1,1)<<init<<make_pair(1,1)<<make_pair(1,1),engine);
		while(!fitdata.AbsoluteOptimalityExitCondition(0.0001))
			fitdata.Iterate(engine);
		PlotFit1D<decltype(fitdata)>()
			.Fit("data_fit_"+suffix,"data_"+suffix,fitdata,0.001)
			<<"set xlabel 'MM, GeV'"<<"set ylabel 'counts'";
		
		double yd=fitdata[0],ym=fitMC[0];
		double dd=fitdata.GetParamParabolicError(0.01,0),dm=fitMC.GetParamParabolicError(0.01,0);
		qBin.y=yd/ym;
		qBin.dy=dd/ym+yd*dm/pow(ym,2);
	}
	PlotHist().Hist("He3eta true events in data",luminosity*=MC_events_count)
		<<"set xlabel 'Q, MeV'"<<"set ylabel 'True events, counts'"<<"set nokey";
	PlotHist lumplot;
	lumplot.Hist("analysed",luminosity/=sigmaHe3eta);
	lumplot.Hist("estimated",luminosity/=PresentRunsAmountRatio("He3"));
	lumplot<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integrated luminocity, 1/nb'";
}
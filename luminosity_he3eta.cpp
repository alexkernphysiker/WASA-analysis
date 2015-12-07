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
#include <str_get.h>
#include <gethist.h>
#include <theory.h>
#include "phys_constants.h"
using namespace std;
using namespace Genetic;
RANDOM engine;
#define Q_range 18.0,28.0
int main(int,char**){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta");
	hist mc_norm(MC,"He3eta",{"Histograms","Reconstruction"},"Reference");
	hist acceptance(MC,"He3eta",{"Histograms","Reconstruction"},"Additional");
	{
		hist mc_filtered1(MC,"He3eta",{"Histograms","Reconstruction"},"Theta_reconstruction_correct");
		hist mc_filtered2(MC,"He3eta",{"Histograms","Reconstruction"},"Reconstructed");
		PlotHist().Hist("All MC events",mc_norm).Hist("FPC",mc_filtered1)
			.Hist("Reconstructed",mc_filtered2).Hist("Preselected",acceptance)
			<<"set yrange [0:]"
			<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	}
	(acceptance/=mc_norm).Cut(Q_range);
	PlotHist().Hist("Acceptance",acceptance)<<"set xlabel 'Q, MeV'"<<"set ylabel 'Acceptance, n.d.'"<<"set nokey";
	hist luminosity=acceptance.CloneEmptyBins();
	LinearInterpolation<double> ChiSqMCPeak,ChiSqData;
	for(auto&qBin:luminosity){
		#define dmmr 0.4,0.6
		#define mmr 0.530,0.553
		int index=int(qBin.x*binning_coefficient);
		hist
			DATAhist(DATA,"He3",{"Histograms","MissingMass"},to_string(index)),
			MC_He3eta(MC,"He3eta",{"Histograms","MissingMass"},to_string(index)),
			MC_He3pi0pi0(MC,"He3pi0pi0",{"Histograms","MissingMass"},to_string(index)),
			MC_He3pi0pi0pi0(MC,"He3pi0pi0pi0",{"Histograms","MissingMass"},to_string(index));
		string suffix=string("Q=")+to_string(qBin.x)+"MeV";
		PlotHist().HistWLine(string("MCHe3eta")+suffix,MC_He3eta.Cut(dmmr))<<"set xlabel 'MM, GeV'"<<"set ylabel 'Counts'";
		PlotHist().HistWLine(string("MCHe3 2pi0")+suffix,MC_He3pi0pi0.Cut(dmmr))<<"set xlabel 'MM, GeV'"<<"set ylabel 'Counts'";
		PlotHist().HistWLine(string("MCHe3 3pi0")+suffix,MC_He3pi0pi0pi0.Cut(dmmr))<<"set xlabel 'MM, GeV'"<<"set ylabel 'Counts'";
		PlotHist().Hist(string("Data")+suffix,DATAhist.Cut(dmmr))<<"set xlabel 'MM, GeV'"<<"set ylabel 'Counts'";
		cout<<qBin.x<<"MeV"<<endl;

		typedef Mul<Par<0>,Func3<Gaussian,Arg<0>,Par<1>,Par<2>>> Peak;
		#define filter make_shared<Above>()<<0<<0<<0
		typedef PolynomFunc<0,3,4> BG;
		const size_t thr=8;
		const double plot_step=0.0001;
		const double accu=0.0000001;

		FitFunction<DifferentialMutations<>,Peak,ChiSquareWithXError> 
			fitMC(make_shared<FitPoints>()<<MC_He3eta.Cut(0.54,0.56));
		fitMC.SetFilter(filter).SetThreadCount(thr)
			.Init(15*Peak::ParamCount,make_shared<GenerateByGauss>()<<make_pair(200,300)<<make_pair(m_eta,0.001)<<make_pair(0.0,0.01),engine);
		while(!fitMC.AbsoluteOptimalityExitCondition(accu))
			fitMC.Iterate(engine);
		cout<<"Peak. ChiSq = "<<fitMC.Optimality()<<endl;
		for(double p:fitMC)cout<<p<<" ";
		cout<<endl;
		PlotFit1D<decltype(fitMC)>().Fit("mc_fit_"+suffix,"mc_"+suffix,fitMC,plot_step)
			<<"set xlabel 'MM, GeV'"<<"set ylabel 'counts'";
		ChiSqMCPeak<<make_pair(qBin.x,fitMC.Optimality());

		FitFunction<DifferentialMutations<>,Add<Peak,BG>,ChiSquareWithXError> 
			fitdata(make_shared<FitPoints>()<<DATAhist.Cut(mmr));
		fitdata.SetFilter(filter).SetThreadCount(thr);
		auto Init=make_shared<GenerateByGauss>()<<make_pair(0,2)<<make_pair(fitMC[1],0.001)<<make_pair(fitMC[2],0.0);
		Init<<make_pair(50000,50000)<<make_pair(-50000,50000)<<make_pair(0,50000);
		while((Init<<make_pair(0,50000))->Count()< BG::ParamCount){}
		fitdata.Init(15*BG::ParamCount,Init,engine);
		while(!fitdata.AbsoluteOptimalityExitCondition(accu))
			fitdata.Iterate(engine);
		cout<<"Data. ChiSq = "<<fitdata.Optimality()<<endl;
		for(double p:fitdata)cout<<p<<" ";
		cout<<endl;
		PlotFit1D<decltype(fitdata)>()
			.Fit("data_fit_"+suffix,"data_"+suffix,fitdata,plot_step)
			.ParamFunc("Peak_"+suffix,Peak(),fitdata,plot_step)
			.ParamFunc("BG_"+suffix,BG(),fitdata,plot_step)
			<<"set xlabel 'MM, GeV'"<<"set ylabel 'counts'";
		ChiSqData<<make_pair(qBin.x,fitdata.Optimality());
		
		double yd=fitdata[0],ym=fitMC[0];
		cout <<"Y = "<<yd<<" ; "<<ym<<endl;
		double dd=fitdata.GetParamParabolicError(0.001,0),dm=fitMC.GetParamParabolicError(0.001,0);
		cout <<"dY = "<<dd<<" ; "<<dm<<endl;
		qBin.y=double(MC_events_count)*yd/ym;
		qBin.dy=double(MC_events_count)*sqrt(pow(dd/ym,2)+pow(dm*yd/pow(ym,2),2));
		#undef mmr
	}
	PlotPoints<double,decltype(ChiSqData)>()
		.LineOnly("Chi^2_data",ChiSqData).LineOnly("Chi^2_MC",ChiSqMCPeak);
	PlotHist().Hist("He3eta true events in data",luminosity)
		<<"set xlabel 'Q, MeV'"<<"set ylabel 'True events, counts'"<<"set nokey";
	PlotHist lastplot;
	lastplot.Hist("analysed",luminosity/=sigmaHe3eta);
	lastplot.Hist("estimated",luminosity/=PresentRunsAmountRatio("He3"))
		<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integrated luminocity, 1/nb'";
}
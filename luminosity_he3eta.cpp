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
using namespace std;
using namespace Genetic;
RANDOM engine;
#define Q_range 15.0,30.0
#define MissingMass_range 0.50,0.57
template<int x_index,int p_index>
class TableData:public virtual IParamFunc{
	LinearInterpolation<double> data;
public:
	TableData(){}
	virtual ~TableData(){}
	virtual double operator()(const ParamSet&X,const ParamSet&P)const override{
		return data(X[x_index])*P[p_index];
	}
	TableData&operator<<(std::pair<double,double>&&point){
		data<<static_cast<std::pair<double,double>&&>(point);
		return *this;
	}
	enum{ParamCount=0,ArgCount=x_index+1};
};
template<int x,int p>
inline std::shared_ptr<TableData<x,p>>&operator<<(
	std::shared_ptr<TableData<x,p>>C,std::pair<double,double>P
){
	C->operator<<(static_cast<std::pair<double,double>&&>(P));
	return C;
}
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
			MC_He3pi0pi0pi0(MC,"He3pi0pi0pi0",{"Histograms","MissingMass"},to_string(index));
		hist DATAhist(DATA,"He3",{"Histograms","MissingMass"},to_string(index));
		DATAhist.Cut(MissingMass_range);
		MC_He3eta.Cut(MissingMass_range);
		MC_He3pi0pi0.Cut(MissingMass_range);
		MC_He3pi0pi0pi0.Cut(MissingMass_range);
		string suffix=string("Q=")+to_string(qBin.x)+"MeV";
		PlotHist().HistWLine(string("MCHe3eta")+suffix,MC_He3eta)
			.HistWLine(string("MCHe3 2pi0")+suffix,MC_He3pi0pi0)
			.HistWLine(string("MCHe3 3pi0")+suffix,MC_He3pi0pi0pi0)
			<<"set xlabel 'MM, GeV'"<<"set ylabel 'Counts'";
		PlotHist().Hist(string("DataHe3")+suffix,DATAhist)<<"set xlabel 'MM, GeV'"<<"set ylabel 'Counts'";
		
		FitFunction<DifferentialMutations<>,
			Add<
				Mul<Func3<Gaussian,Arg<0>,Par<1>,Par<2>>,Par<0>>,
				Add<TableData<0,3>,TableData<0,4>>
			>
			,ChiSquareWithXError> fit(make_shared<FitPoints>()<<DATAhist);
		for(auto&P:MC_He3pi0pi0)
			dynamic_pointer_cast<TableData<0,3>>(fit.Func())<<make_pair(P.x,P.y);
		for(auto&P:MC_He3pi0pi0pi0)
			dynamic_pointer_cast<TableData<0,4>>(fit.Func())<<make_pair(P.x,P.y);
	}
	const double MC_events_count=5000000;
	PlotHist().Hist("He3eta true events in data",luminosity*=MC_events_count)
		<<"set xlabel 'Q, MeV'"<<"set ylabel 'True events, counts'"<<"set nokey";
	PlotHist lumplot;
	lumplot.Hist("analysed",luminosity/=sigmaHe3eta);
	lumplot.Hist("estimated",luminosity/=PresentRunsAmountRatio("He3"));
	lumplot<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integrated luminocity, 1/nb'";
}
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

template<int p,class FUNC>
class TableData:public virtual FUNC{
	LinearInterpolation<double> data;
public:
	TableData(){}
	virtual ~TableData(){}
	virtual double operator()(const ParamSet&X,const ParamSet&P)const override{
		return data(FUNC::operator()(X,P))*P[p];
	}
	TableData&operator<<(std::pair<double,double>&&point){
		data<<static_cast<std::pair<double,double>&&>(point);
		return *this;
	}
	size_t size()const{return data.size();}
	std::pair<double,double>&operator[](size_t i)const{return data[i];}
	enum{ParamCount=max2<FUNC::ParamCount,p+1>::val,ArgCount=FUNC::ArgCount};
};
template<int p,class FUNC>
inline std::shared_ptr<TableData<p,FUNC>> operator<<(
	std::shared_ptr<TableData<p,FUNC>>C,std::pair<double,double>P
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
#define init make_shared<GenerateByGauss>()<<make_pair(100,100)<<make_pair(m_eta,0.01)<<make_pair(0.01,0.01)
#define filter make_shared<Above>()<<0<<0<<0
		fitMC.SetFilter(filter).Init(pop_size,init,engine);
		while(!fitMC.RelativeParametersDispersionExitCondition(0.0001))
			fitMC.Iterate(engine);
		PlotFit1D<decltype(fitMC)>()
			.Fit("mc_fit_"+suffix,"mc_"+suffix,fitMC,0.00001)
			<<"set xlabel 'MM, GeV'"<<"set ylabel 'counts'";
	}
	PlotHist().Hist("He3eta true events in data",luminosity*=MC_events_count)
		<<"set xlabel 'Q, MeV'"<<"set ylabel 'True events, counts'"<<"set nokey";
	PlotHist lumplot;
	lumplot.Hist("analysed",luminosity/=sigmaHe3eta);
	lumplot.Hist("estimated",luminosity/=PresentRunsAmountRatio("He3"));
	lumplot<<"set xlabel 'Q, MeV'"<<"set ylabel 'Integrated luminocity, 1/nb'";
}
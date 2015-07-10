// this file is distributed under 
// GPL v 3.0 license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>
#include <fit.h>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
#include <phys_constants.h>
#include "gethist.h"
using namespace std;
using namespace Genetic;
typedef LinearInterpolation<double> FuncTbl;
typedef PlotPoints<double,FuncTbl> PlotTbl;
#define ParR(x) static_cast<ParamSet&&>(x)
int main(int,char**){
#include "env.cc"
	Plotter::Instance().SetOutput(outpath+"/He3Eta");
	vector<string> Kinematics={"Histograms","Kinematics"};
	vector<string> EventsCount={"Histograms","EventsCount"};
	pair<FuncTbl,FuncTbl> acceptance;
	vector<pair<double,pair<FuncTbl,FuncTbl>>> missingmass_mc_normed;
	{
		string MCFile=inputpath+"/MCHe3eta_gg_.root";
		pair<FuncTbl,FuncTbl> TotalEvents;{
			hist normhist(MCFile,EventsCount,"AllEventsOnPBeam");
			for(auto &p:normhist){
				TotalEvents.first<<make_pair(p.x,p.y);
				TotalEvents.second<<make_pair(p.x,p.dy);
			}
		}
		hist registered(MCFile,EventsCount,"FilteredEventsOnPBeam");
		PlotTbl missingmass_spectra_plot;
		for(auto &histpoint:registered){
			double norm=TotalEvents.first(histpoint.x);
			double dnorm=TotalEvents.second(histpoint.x);
			if(norm>48000){
				hist subhist(MCFile,Kinematics,string("MissingMass")+to_string(int(histpoint.x*1000)));
				pair<FuncTbl,FuncTbl> spectrum;
				for(auto &subhistpoint:subhist){
					spectrum.first<<make_pair(subhistpoint.x,subhistpoint.y/norm);
					spectrum.second<<make_pair(subhistpoint.x,(dnorm*subhistpoint.y/pow(norm,2))+(subhistpoint.dy/norm));
				}
				missingmass_spectra_plot.WithErrorOnX(
					"MC missing mass P="+to_string(histpoint.x),static_cast<FuncTbl&&>(spectrum.first),
					[&spectrum](double x){return spectrum.second(x);}
				);
				double acc=histpoint.y/norm;
				if(acc>0.1){
					missingmass_mc_normed.push_back(make_pair(histpoint.x,spectrum));
					acceptance.first<<make_pair(histpoint.x,acc);
					acceptance.second<<make_pair(histpoint.x,(dnorm*histpoint.y/pow(norm,2))+(histpoint.dy/norm));
				}
			}
		}
		PlotTbl().WithErrorOnX(
			"Acceptance",static_cast<FuncTbl&&>(acceptance.first),
			[&acceptance](double x){return acceptance.second(x);}
		);
	}
	vector<pair<double,shared_ptr<FitPoints>>> missing_mass_data;{
		vector<pair<double,hist>> missing_mass_spectra;
		for(int run_no=46200,inited=0;run_no<=46250;run_no++){
			string rootfile=inputpath+"/DataHe3eta_gg__run_"+to_string(run_no)+".root";
			hist alleventshist(rootfile,EventsCount,"AllEventsOnPBeam");
			if(alleventshist.count()>0){
				printf("File %s found...\n",rootfile.c_str());
				for(auto &histparam:missingmass_mc_normed){
					double p_beam=histparam.first;
					hist subhist(rootfile,Kinematics,string("MissingMass")+to_string(int(p_beam*1000)));
					if(0==inited)
						missing_mass_spectra.push_back(make_pair(p_beam,subhist));
					else{
						for(pair<double,hist>&P:missing_mass_spectra)
							if(P.first==p_beam)
								P.second+=subhist;
					}
				}
				inited=1;
			}
		}
		for(auto &P:missing_mass_spectra){
			double p_beam=P.first;
			auto points=make_shared<FitPoints>();
			for(hist::point p: P.second){
				FitPoints::Point point;
				point.X<<p.x;
				point.WX<<p.dx;
				point.y=p.y;
				point.wy=p.dy;
				points<<point;
			}
			missing_mass_data.push_back(make_pair(p_beam,points));
		}
	}
	pair<FuncTbl,FuncTbl> Events_Count;
	for(auto &spectrum: missing_mass_data){
		double p_beam=spectrum.first;
		shared_ptr<FitPoints> points=spectrum.second;
		string suffix="P="+to_string(p_beam);
		
		printf("Missing mass spectrum %s...\n",suffix.c_str());
		PolynomFunc<0,1,4> BackGround;
		pair<FuncTbl,FuncTbl> ForeGround;
		pair<paramFunc,paramFunc> fitfuncs=make_pair(
			[&ForeGround,&BackGround](ParamSet&&X,ParamSet&&P){
				return (P[0]*ForeGround.first(X[0]))+BackGround(ParR(X),ParR(P));
			},
			[&ForeGround](ParamSet&&X,ParamSet&&P){
				return P[0]*ForeGround.second(X[0]);
			}
		);
		for(auto &norm_hist:missingmass_mc_normed)
			if(norm_hist.first==p_beam)
				ForeGround=norm_hist.second;
		printf("norm. cnt= %i and %i \n",ForeGround.first.size(),ForeGround.second.size());
		
		PlotPoints1D().Points("Missing mass DATA "+suffix,points);
		
		auto pointstofit=SelectFitPoints(points,make_shared<Filter>([](ParamSet&&X){return abs((m_eta-X[0])*1000)<50;}));
		printf("Preparing fit...points count = %i of %i \n",pointstofit->count(),points->count());
		FitFunctionWithError<Crossing<DifferentialMutations<>>,ChiSquareWithXError> fit(
			pointstofit,fitfuncs.first,fitfuncs.second
		);
		fit.SetFilter([](ParamSet&&P){return P[0]>0;});
		fit.SetCrossingProbability(0.01);
		auto init=make_shared<GenerateByGauss>()<<make_pair(1000,1000)<<make_pair(0,5000)<<make_pair(0,1000);
		while(init->Count()<BackGround.ParamCount)
			init<<make_pair(0,0.01);
		fit.Init(BackGround.ParamCount*15,init);
		
		printf("Fitting...\n");
		while(!fit.AbsoluteOptimalityExitCondition(0.000001)){
			fit.Iterate();
			printf("%i iterations. %f<=chi^2<=%f   \r",fit.iteration_count(),fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
		}
		printf("\n");
		Events_Count.first<<make_pair(p_beam,fit[0]);
		Events_Count.second<<make_pair(p_beam,fit.GetParamParabolicError(0.01,0));
		PlotFit1D<decltype(fit)>().Points("Data",pointstofit).Fit("Fit",fit,0.01)
			.ParamFunc("BG",static_cast<IParamFunc&&>(BackGround),fit,0.01)
			.ParamFunc("FG",[&ForeGround](ParamSet&&X,ParamSet&&P){return P[0]*ForeGround.first(X[0]);},fit,0.01);
	}
	PlotTbl().WithErrorOnX("Events count",static_cast<FuncTbl&&>(Events_Count.first),[&Events_Count](double x){return Events_Count.second(x);});
}
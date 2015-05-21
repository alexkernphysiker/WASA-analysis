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
#include "gnuplot.h"
using namespace std;
using namespace Genetic;
typedef LinearInterpolation<double> FuncTbl;
#define ParR(x) static_cast<ParamSet&&>(x)
int main(int,char**){
#include "env.cc"
	SetPlotOutput(outpath+"/He3Eta");
	vector<string> Kinematics;
	Kinematics.push_back("Histograms");
	Kinematics.push_back("Kinematics");
	vector<string> EventsCount;
	EventsCount.push_back("Histograms");
	EventsCount.push_back("EventsCount");
	pair<FuncTbl,FuncTbl> acceptance;
	vector<pair<double,pair<FuncTbl,FuncTbl>>> missingmass_mc_normed;
	{
		string MCFile=inputpath+"/MCHe3Eta.root";
		pair<FuncTbl,FuncTbl> TotalEvents;{
			hist normhist(MCFile,EventsCount,"AllEventsOnPBeam");
			for(auto &p:normhist){
				TotalEvents.first<<make_pair(p.x,p.y);
				TotalEvents.second<<make_pair(p.x,p.dy);
			}
		}
		hist registered(MCFile,EventsCount,"FilteredEventsOnPBeam");
		Plot missingmass_spectra_plot;
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
				missingmass_spectra_plot.Points(
					"MC missing mass P="+to_string(histpoint.x),
					spectrum.first,[&spectrum](double x){return spectrum.second(x);}
				);
				double acc=histpoint.y/norm;
				if(acc>0.1){
					missingmass_mc_normed.push_back(make_pair(histpoint.x,spectrum));
					acceptance.first<<make_pair(histpoint.x,acc);
					acceptance.second<<make_pair(histpoint.x,(dnorm*histpoint.y/pow(norm,2))+(histpoint.dy/norm));
				}
			}
		}
		Plot acceptance_plot;
		acceptance_plot.Points("Acceptance",acceptance.first,[&acceptance](double x){return acceptance.second(x);});
	}
	vector<pair<double,shared_ptr<FitPoints>>> missing_mass_data;{
		vector<pair<double,hist>> missing_mass_spectra;
		for(int run_no=46200,inited=0;run_no<=46250;run_no++){
			string rootfile=inputpath+"/DataHe3Eta_run_"+to_string(run_no)+".root";
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
		for(auto &norm_hist:missingmass_mc_normed)
			if(norm_hist.first==p_beam)
				ForeGround=norm_hist.second;
		printf("norm. cnt= %i and %i \n",ForeGround.first.size(),ForeGround.second.size());
		
		Plot fitplot;
		fitplot.Points("Missing mass DATA "+suffix,points);
		
		auto pointstofit=SelectFitPoints(points,make_shared<Filter>([](ParamSet&&X){return abs((m_eta-X[0])*1000)<50;}));
		printf("Preparing fit...points count = %i of %i \n",pointstofit->count(),points->count());
		FitFunctionWithError<Crossing<DifferentialMutations<>>,ChiSquareWithXError> fit(
			pointstofit,
			[&ForeGround,&BackGround](ParamSet&&X,ParamSet&&P){
				return (P[0]*ForeGround.first(X[0]))+BackGround(ParR(X),ParR(P));
			},
			[&ForeGround](ParamSet&&X,ParamSet&&P){
				return P[0]*ForeGround.second(X[0]);
			}
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
		fitplot.Function("fit "+suffix,
			[&fit](double x){return fit(ParamSet(x));},
			pointstofit->operator[](0).X[0],pointstofit->operator[](pointstofit->count()-1).X[0],0.01
		);
		fitplot.Function("Background "+suffix,
			[&fit,&BackGround](double x){return BackGround(ParamSet(x),fit.Parameters());},
			pointstofit->operator[](0).X[0],pointstofit->operator[](pointstofit->count()-1).X[0],0.01
		);
		fitplot.Function("Foreground "+suffix,
			[&fit,&ForeGround](double x){return fit[0]*ForeGround.first(x);},
			pointstofit->operator[](0).X[0],pointstofit->operator[](pointstofit->count()-1).X[0],0.01
		);
	}
	Plot eventsplot;
	eventsplot.Points("Events count",Events_Count.first,[&Events_Count](double x){return Events_Count.second(x);});
}
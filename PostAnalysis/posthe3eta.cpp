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
			for(hist::point p:normhist){
				TotalEvents.first<<make_pair(p.x,p.y);
				TotalEvents.second<<make_pair(p.x,p.dy);
			}
		}
		hist registered(MCFile,EventsCount,"FilteredEventsOnPBeam");
		Plot missingmass_spectra_plot;
		for(hist::point histpoint:registered){
			double norm=TotalEvents.first(histpoint.x);
			double dnorm=TotalEvents.second(histpoint.x);
			hist subhist(MCFile,Kinematics,string("MissingMass")+to_string(int(histpoint.x*1000)));
			pair<FuncTbl,FuncTbl> spectrum;
			for(hist::point subhistpoint:subhist){
				spectrum.first<<make_pair(subhistpoint.x,subhistpoint.y/norm);
				spectrum.second<<make_pair(subhistpoint.x,(dnorm*subhistpoint.y/pow(norm,2))+(subhistpoint.dy/norm));
			}
			missingmass_spectra_plot.Points(
				"MC missing mass P="+to_string(histpoint.x),
				spectrum.first,[&spectrum](double x){return spectrum.second(x);});
			missingmass_mc_normed.push_back(make_pair(histpoint.x,spectrum));
			acceptance.first<<make_pair(histpoint.x,histpoint.y/norm);
			acceptance.second<<make_pair(histpoint.x,(dnorm*histpoint.y/pow(norm,2))+(histpoint.dy/norm));
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
				for(auto histparam:missingmass_mc_normed){
					double p_beam=histparam.first;
					hist subhist(rootfile,Kinematics,string("MissingMass")+to_string(int(p_beam*1000)));
					if(0==inited)
						missing_mass_spectra.push_back(make_pair(p_beam,subhist));
					else{
						for(auto P:missing_mass_spectra)
							if(P.first==p_beam)
								P.second+=subhist;
					}
				}
				inited=1;
			}
		}
		for(auto P:missing_mass_spectra){
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
	for(auto spectrum: missing_mass_data){
		double p_beam=spectrum.first;
		shared_ptr<FitPoints> points=spectrum.second;
		Plot fitplot;
		string suffix=to_string(int(p_beam*1000));
		fitplot.Points("Missing mass DATA "+suffix,points);
		PolynomFunc<0,3,6> BackGround;
		pair<FuncTbl,FuncTbl> *ForeGround=nullptr;
		for(auto norm_hist:missingmass_mc_normed)
			if(norm_hist.first==p_beam)
				ForeGround=&norm_hist.second;
		function<double(ParamSet&,ParamSet&)> F;
		function<double(ParamSet&,ParamSet&)> E;
		Fit2<DifferentialMutations<>,ChiSquareWithXError> fit(
			points,
			[ForeGround,&BackGround](ParamSet&X,ParamSet&P){
				return (P[0]*ForeGround->first(P[1]*X[0]+P[2]))+BackGround(X,P);
			},
			[ForeGround](ParamSet&X,ParamSet&P){
				return P[0]*ForeGround->second(P[1]*X[0]+P[2]);
			}
		);
	}
}
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>
#include <paramfunc.h>
#include <fitpoints.h>
#include <filter.h>
#include <initialconditions.h>
#include <genetic.h>
#include <phys_constants.h>
#include "gethist.h"
#include "gnuplot.h"
using namespace std;
using namespace Fit;
typedef LinearInterpolation<double> FuncTbl;
typedef Mul<Par<0>,Func3<Gaussian,Arg<0>,Par<1>,Par<2>>> Foreground;
int main(int,char**){
#include "env.cc"
	SetPlotOutput(outpath+"/He3Eta");
	string MCFile=inputpath+"/MCHe3Eta.root";
	string DataFile=inputpath+"/DataHe3Eta.root";
	vector<string> Kinematics;
	Kinematics.push_back("Histograms");
	Kinematics.push_back("Kinematics");
	vector<string> EventsCount;
	EventsCount.push_back("Histograms");
	EventsCount.push_back("EventsCount");
	FuncTbl acceptance,dacceptance;
	vector<pair<double,pair<FuncTbl,FuncTbl>>> missingmass_mc_normed;
	{
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
				double dy=subhistpoint.dy;
				if(dy<1)dy=1;
				spectrum.second<<make_pair(subhistpoint.x,(dnorm*subhistpoint.y/pow(norm,2))+(dy/norm));
			}
			missingmass_spectra_plot.Points(
				"MC missing mass P="+to_string(histpoint.x),
				spectrum.first,[&spectrum](double x){return spectrum.second(x);});
			missingmass_mc_normed.push_back(make_pair(histpoint.x,spectrum));
			acceptance<<make_pair(histpoint.x,histpoint.y/norm);
			dacceptance<<make_pair(histpoint.x,(dnorm*histpoint.y/pow(norm,2))+(histpoint.dy/norm));
		}
		Plot acceptance_plot;
		acceptance_plot.Points("Acceptance",acceptance,[&dacceptance](double x){return dacceptance(x);});
	}
}
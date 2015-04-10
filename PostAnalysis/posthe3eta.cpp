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
typedef Mul<Par<0>,Func3<Gaussian,Arg<0>,Par<1>,Par<2>>> Foreground;
int main(int,char**){
#include "env.cc"
	SetPlotOutput(outpath+"/He3Eta");
	string MCFile=inputpath+"/MCHe3Eta.root";
	vector<string> Kinematics;
	Kinematics.push_back("Histograms");
	Kinematics.push_back("Kinematics");
	vector<string> EventsCount;
	EventsCount.push_back("Histograms");
	EventsCount.push_back("EventsCount");
	LinearInterpolation<double> acceptance,dacceptance;
	vector<pair<double,LinearInterpolation<double>>> missingmass_mc_normed;
	{
		LinearInterpolation<double> mc_norm,mc_dnorm;
		{
			hist normhist(MCFile,EventsCount,"AllEventsOnPBeam");
			for(hist::point p:normhist){
				mc_norm<<make_pair(p.x,p.y);
				mc_dnorm<<make_pair(p.x,p.dy);
			}
		}
		hist registered(MCFile,EventsCount,"FilteredEventsOnPBeam");
		for(hist::point histpoint:registered){
			double norm=mc_norm(histpoint.x);
			double dnorm=mc_dnorm(histpoint.x);
			if(mc_norm(histpoint.x)>50000){
				string subhistname=string("MissingMass")+to_string(int(histpoint.x*1000));
				hist subhist(MCFile,Kinematics,subhistname);
				LinearInterpolation<double> points;
				for(hist::point subhistpoint:subhist)
					points<<make_pair(subhistpoint.x,subhistpoint.y/norm);
				missingmass_mc_normed.push_back(make_pair(histpoint.x,points));
				acceptance<<make_pair(histpoint.x,histpoint.y/norm);
				dacceptance<<make_pair(histpoint.x,(dnorm*histpoint.y/pow(norm,2))+(histpoint.dy/norm));
			}
		}
		Plot fig1;
		fig1.Points("Acceptance",acceptance,[&dacceptance](double x){return dacceptance(x);});
	}
}
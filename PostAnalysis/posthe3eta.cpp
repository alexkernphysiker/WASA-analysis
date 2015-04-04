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
int main(int,char**) {
#include "env.cc"
	string MCFile=inputpath+"/MCHe3Eta.root";
	vector<string> Kinematics;
	Kinematics.push_back("Histograms");
	Kinematics.push_back("Kinematics");
	vector<string> Reconstruction;
	Reconstruction.push_back("Histograms");
	Reconstruction.push_back("Reconstruction");
	vector<string> EventsCount;
	EventsCount.push_back("Histograms");
	EventsCount.push_back("OutputEventsCount");
	Plot fig1(outpath);
	LinearInterpolation<double> norm,dnorm;{
		hist normhist(MCFile,Reconstruction,"P_beam");
		for(point p:normhist){
			norm<<make_pair(p.x,p.y);
			dnorm<<make_pair(p.x,p.dy);
		}
	}
	LinearInterpolation<double> acceptance,dacceptance;{
		hist registered(MCFile,EventsCount,"DependenceOnBeam");
		for(point p:registered)
			if(norm(p.x)>0){
				acceptance<<make_pair(p.x,p.y/norm(p.x));
				dacceptance<<make_pair(p.x,(dnorm(p.x)*p.y/pow(norm(p.x),2))+(p.dy/norm(p.x)));
				printf("%f\n",p.x);
			}
	}
	printf("plotting\n");
	fig1.Points("Acceptance",acceptance,[&dacceptance](double x){return dacceptance(x);});
	fig1.Out("figure1",true);
}
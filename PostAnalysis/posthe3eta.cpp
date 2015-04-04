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
	hist MissingMass_all(MCFile,Kinematics,"MissingMass_all");
	fig1.Hist("obtained_data",MissingMass_all);
	fig1.Out("fig1");
}
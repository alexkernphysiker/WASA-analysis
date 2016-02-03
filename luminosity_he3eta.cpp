// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <math_h/functions.h>
#include <Genetic/fit.h>
#include <Genetic/paramfunc.h>
#include <Genetic/filter.h>
#include <Genetic/initialconditions.h>
#include <str_get.h>
#include <gethist.h>
#include <theory.h>
#include "phys_constants.h"
#include "kinematics.h"
using namespace std;
using namespace ROOT_data;
using namespace GnuplotWrap;
#define Q_range 18.0,30.0
int main(int,char**){
	Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta forward");
	vector<string> histpath_forward={"Histograms","He3Forward_Reonstruction"};
	hist mc_norm(MC,"He3eta",histpath_forward,"0-Reference");
	double MC_events_count=0;
	for(const hist::point&P: mc_norm)MC_events_count+=P.y;
	cout<<"Montecarlo evens count of "<<MC_events_count<<" detected."<<endl;
	{
		hist n1(MC,"He3eta",histpath_forward,"1-AllTracks");
		hist n2(MC,"He3eta",histpath_forward,"2-FPC");
		hist n3(MC,"He3eta",histpath_forward,"3-AllCuts");
		hist n4(MC,"He3eta",histpath_forward,"4-Reconstructed");
		PlotHist().Hist("All MC events",mc_norm).Hist("Forward Tracks",n1).Hist("FPC",n2).Hist("E_{dep} cuts",n3)
			<<"set yrange [0:]"<<"set xlabel 'Q, MeV'"<<"set ylabel 'Events count'";
	}
	//ToDo: perform analysis due to new preselection data format
}
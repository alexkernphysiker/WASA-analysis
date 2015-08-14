// this file is distributed under 
// GPL v 3.0 license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <fit.h>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
#include "gethist.h"
#include "fit_hist2hist.h"
using namespace std;
using namespace Genetic;
RANDOM engine;
int main(int,char**){
#include "env.cc"
	Plotter::Instance().SetOutput(outpath,"he3eta");
	vector<string> kin_path={"Histograms","Kinematics"};
	vector<string> norm_path={"Histograms","EventsCount"};
	
	auto mc_norm=make_shared<hist>(false,"He3eta",static_right(norm_path),"AllEventsOnPBeam");
	auto mc_filtered=make_shared<hist>(false,"He3eta",static_right(norm_path),"FilteredEventsOnPBeam");
	auto mc_reconstructed=make_shared<hist>(false,"He3eta",static_right(norm_path),"ReconstructedEventsOnPBeam");
	PlotHist().Hist("All MC events",mc_norm).Hist("Preselected",mc_filtered).Hist("Reconstructed",mc_reconstructed);
	
	auto acceptance=make_shared<hist>(*mc_reconstructed);acceptance->operator/=(*mc_norm);
	PlotHist().Hist("Acceptance",acceptance);
	
	for(hist::point BeamMomentaBin:*acceptance){
		
	}
}
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
const size_t bg_power=2;
int main(int,char**){
#include "env.cc"
	Plotter::Instance().SetOutput(outpath,"he3eta");
	vector<string> kin_path={"Histograms","Kinematics"};
	auto missingmass_mc=make_shared<hist>(false,"He3eta_gg_",static_right(kin_path),"MissingMass1633");
	auto missingmass_data=make_shared<hist>(true,"He3eta_gg_",static_right(kin_path),"MissingMass1633");
	FitHist<DifferentialMutations<>> fit(missingmass_data,missingmass_mc,[](double x,ParamSet&&P){
		return Polynom(x,P,bg_power,1);
	});
	PlotHist().Hist("MC",missingmass_mc);
	PlotHist().Hist("Data",missingmass_data);
}
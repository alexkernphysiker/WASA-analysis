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
Hist2Hist::bg_func bg_polynom=[](double x,ParamSet&&P){return Polynom(x,P,bg_power,1);};
int main(int,char**){
#include "env.cc"
	Plotter::Instance().SetOutput(outpath,"he3eta");
	vector<string> kin_path={"Histograms","Kinematics"};
	auto missingmass_mc=make_shared<hist>(false,"He3eta_gg_",static_right(kin_path),"MissingMass1633");
	auto missingmass_data=make_shared<hist>(true,"He3eta_gg_",static_right(kin_path),"MissingMass1633");
	FitHist<DifferentialMutations<>> fit(missingmass_data,missingmass_mc,bg_polynom);
	auto init=make_shared<GenerateByGauss>()<<make_pair(1,1);
	for(size_t i=0;i<=bg_power;i++)init<<make_pair(0,1);
	fit.Init((bg_power+2)*20,init);
	
	PlotHist().Hist("MC",missingmass_mc);
	PlotHist().Hist("Data",missingmass_data);
}
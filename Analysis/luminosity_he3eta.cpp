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
using namespace std;
using namespace Genetic;
int main(int,char**){
#include "env.cc"
	Plotter::Instance().SetOutput(outpath,"he3eta");
	hist missingmass_mc(false,"He3eta_gg_",{"Histograms","Kinematics"},"MissingMass1633");
	hist missingmass_data(true,"He3eta_gg_",{"Histograms","Kinematics"},"MissingMass1633");
	PlotHist().Hist("MC",static_right(missingmass_mc));
	PlotHist().Hist("Data",static_right(missingmass_data));
}
// this file is distributed under 
// GPL v 3.0 license
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
using namespace std;
using namespace Genetic;
typedef LinearInterpolation<double> FuncTbl;
typedef PlotPoints<double,FuncTbl> PlotTbl;
#define ParR(x) static_cast<ParamSet&&>(x)
int main(int,char**){
#include "env.cc"
	string simulation=inputpath+"../Reconstruction";
	Plotter::Instance().SetOutput(outpath+"../Reconstruction");
}
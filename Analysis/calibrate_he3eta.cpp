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
#include "read_simulation.h"
using namespace std;
using namespace Genetic;
typedef LinearInterpolation<double> FuncTbl;
typedef PlotPoints<double,FuncTbl> PlotTbl;
#define ParR(x) static_cast<ParamSet&&>(x)
int main(int,char**){
#include "env.cc"
	string simulation=inputpath+"../Reconstruction/";
	Plotter::Instance().SetOutput(outpath+"../Reconstruction/");
	{//Energy
		printf("\nEnergy\n");
		auto data=ReadWeightedFrom2D(
			inputpath+"He3.E.simulation.txt",
			0.08,0.024,50,
			0.22,0.4,50,
			[](double&,double&){return true;}
		);
		FitFunction<DifferentialMutations<>,PolynomFunc<0,0,2>,SumWeightedSquareDiff> fit(data);
		fit.Init(50,make_shared<GenerateByGauss>()<<make_pair(0,0.1)<<make_pair(1,0.1)<<make_pair(0,0.1));
		while(!fit.RelativeOptimalityExitCondition(0.001)){
			fit.Iterate();
			printf("%i iterations. %f<=S<=%f\r",fit.iteration_count(),fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
		}
		PlotFit1D<decltype(fit)>().Points("points-E",data).Fit("He3.E.calibration.txt",fit,0.001);
	}
	{//phi
		printf("\nPhi\n");
		auto data=ReadWeightedFrom2D(
			inputpath+"He3.phi.simulation.txt",
			0.0,6.3,50,
			0.0,6.3,50,
			[](double&x,double&y){return pow(x-y,2)<0.5;}
		);
		FitFunction<DifferentialMutations<>,PolynomFunc<0,0,2>,SumWeightedSquareDiff> fit(data);
		fit.Init(50,make_shared<GenerateByGauss>()<<make_pair(0,0.1)<<make_pair(1,0.1)<<make_pair(0,0.1));
		while(!fit.RelativeOptimalityExitCondition(0.001)){
			fit.Iterate();
			printf("%i iterations. %f<=S<=%f\r",fit.iteration_count(),fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
		}
		PlotFit1D<decltype(fit)>().Points("points-E",data).Fit("He3.phi.calibration.txt",fit,0.001);
	}
}
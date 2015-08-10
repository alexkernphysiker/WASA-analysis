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
RANDOM engine;
int main(int,char**){
#include "env.cc"
	Plotter::Instance().SetOutput(outpath,"he3eta");
	vector<string> kin_path={"Histograms","Kinematics"};
	auto missingmass_mc=make_shared<hist>(false,"He3eta",static_right(kin_path),"MissingMass1633");
	auto missingmass_data=make_shared<hist>(true,"He3eta",static_right(kin_path),"MissingMass1633");
	FitHist<DifferentialMutations<>> fit(missingmass_data,missingmass_mc,[](double x,ParamSet&&P){return 0.0;});
	auto init=make_shared<GenerateByGauss>()<<make_pair(0.5,0.5);
	fit.SetFilter([](ParamSet&&P){return P[0]>0;});
	printf("%i par;\n",init->Count());
	fit.Init(init->Count()*30,init,engine);
	while(!fit.RelativeOptimalityExitCondition(0.000001)){
		fit.Iterate(engine);
		printf("%i;%f<=S<=%f        \r",fit.iteration_count(),fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
	}
	printf("\n");
	for(double p:fit)printf("%f; ",p);
	printf("\n");
	PlotHist().Hist("MC",missingmass_mc);
	LinearInterpolation<double> f=fit.GetFitFunction(),fg=fit.GetForeground(),bg=fit.GetBackground();
	PlotHist().Hist("Data",missingmass_data).
		Line("fit",[&f](double x){return f(x);},0.5,0.6,0.001).
		Line("fg",[&fg](double x){return fg(x);},0.5,0.6,0.001).
		Line("bg",[&bg](double x){return bg(x);},0.5,0.6,0.001);
}
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
	auto missingmass_he3eta=make_shared<hist>(false,"He3eta",static_right(kin_path),"MissingMass1633");
	PlotHist().Hist("MC he3 eta",missingmass_he3eta);
	auto missingmass_he3pi0pi0=make_shared<hist>(false,"He3pi0pi0",static_right(kin_path),"MissingMass1633");
	PlotHist().Hist("MC he3 2pi0",missingmass_he3pi0pi0);
	auto missingmass_he3pi0pi0pi0=make_shared<hist>(false,"He3pi0pi0pi0",static_right(kin_path),"MissingMass1633");
	PlotHist().Hist("MC he3 3pi0",missingmass_he3pi0pi0pi0);
	auto missingmass_data=make_shared<hist>(true,"He3eta",static_right(kin_path),"MissingMass1633");
	HistFitWithShift<DifferentialMutations<>,HistToHists,shared_ptr<hist>,vector<shared_ptr<hist>>&&> 
		fit(missingmass_data,{missingmass_he3eta});
	printf("initing...\n");
	auto init=make_shared<GenerateByGauss>()<<make_pair(0,0)<<make_pair(0.5,0.5);
	fit.SetFilter([](ParamSet&&P){
		for(size_t i=1;i<P.Count();i++)if(P[i]<0)return false;
		return true;
	});
	printf("%i par;\n",init->Count());
	fit.Init(init->Count()*30,init,engine);
	printf("fitting...\n");
	while(!fit.RelativeOptimalityExitCondition(0.0000001)){
		fit.Iterate(engine);
		printf("%i;%f<=S<=%f        \r",fit.iteration_count(),fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
	}
	printf("\n");
	for(double p:fit)printf("%f; ",p);
	printf("\n");
	PlotHist().Hist("Data",missingmass_data).Hist("eta",fit.GetOpt().GetSubHist(0,fit.Parameters()));
}
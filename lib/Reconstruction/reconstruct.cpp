// this file is distributed under 
// GPL v 3.0 license
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <memory>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
#include <str_get.h>
#include "read_simulation.h"
#include "reconstruct.h"
using namespace std;
using namespace Genetic;
string SimulationDataPath(){
	static string str="";
	if(str=="")str=ENV(PRESEL_DATA)+string("/../Reconstruction/");
	return str;
}
void ProcessReconstruction(string&&name,double from,double to, size_t bins,function<bool(double&,double&)> cut,RANDOM&engine){
	printf("\nProcessing \"%s\"...\n",name.c_str());
	auto data=ReadWeightedFrom2D(SimulationDataPath()+name+".simulation.txt",from,to,bins,from,to,bins,cut);
	printf("%i points with weights\n",data->count());
	typedef PolynomFunc<0,0,3> fnc;
	FitFunction<DifferentialMutations<>,fnc,SumWeightedSquareDiff> fit(data);
	auto init=make_shared<GenerateByGauss>()<<make_pair(0,0.1)<<make_pair(1,0.1);
	while(init->Count()<fnc::ParamCount)init<<make_pair(0,0.01);
	fit.Init(20*fnc::ParamCount,init,engine);
	printf("Fitting...\n");
	while(!fit.RelativeOptimalityExitCondition(0.000001)){
		fit.Iterate(engine);
		printf("%i iterations. %f<=S<=%f               \r",fit.iteration_count(),fit.Optimality(),fit.Optimality(fit.PopulationSize()-1));
	}
	printf("\n");
	Plot<double>().File(name+".simulation.txt",name+" Points","using 1:2");
	Plot<double>()
	.File(name+".simulation.txt.cut.txt",name+" Points","using 1:2")
	.Line(name+".calibration",[&fit](double x){return fit(ParamSet(x));},from,to,0.001);
	ofstream file;
	file.open((name+".params.txt").c_str());
	for(double p:fit)file<<p<<" ";
	file.close();
}
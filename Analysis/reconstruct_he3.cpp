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
#include "read_simulation.h"
using namespace std;
using namespace Genetic;
RANDOM engine;
int main(int,char**){
#include "env.cc"
	string simulation=inputpath+"/../Reconstruction/";
	Plotter::Instance().SetOutput(outpath+"/../Reconstruction/","he3reconstruction");
	ofstream file;
	file.open(simulation+"he3eta.parameters.txt");
	if(file.is_open()){
		auto processfile=[&simulation,&file](string&&name,double from,double to, size_t bins,function<bool(double&,double&)> cut){
			printf("\nProcessing \"%s\"...\n",name.c_str());
			auto data=ReadWeightedFrom2D(simulation+name+".simulation.txt",from,to,bins,from,to,bins,cut);
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
			file<<name.c_str()<<"\n";
			for(double p:fit)file<<p<<" ";
			file<<"\n";
		};
		
		processfile("He3.E.FTH1",0.00,0.20,40,[](double&,double&){return true;});
		processfile("He3.E.FRH1",0.00,0.60,60,[](double&,double&){return true;});
		processfile("He3.E.FRH2",0.25,0.60,60,[](double&x,double&y){return (x>0.25)&&(y>0.45);});
		
		processfile("He3.th",0,0.2,60,[](double&x,double&y){return pow(x-y-0.01,2)<0.0003;});
		processfile("He3.phi",0,6.3,100,[](double&x,double&y){return pow(x-y,2)<0.5;});
		file.close();
	}
}
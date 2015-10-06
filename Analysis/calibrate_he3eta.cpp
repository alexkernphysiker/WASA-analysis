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
		{//Energy
			printf("\nEnergy1...\n");
			auto data=ReadWeightedFrom2D(
				simulation+"He3.E.FRH1.simulation.txt",
				0.0,0.5,100,
				0.0,0.5,100,
				[](double&,double&){return true;}
			);
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
			Plot<double>().Line("He3.E.FRH1.calibration",[&fit](double x){return fit(ParamSet(x));},0.0,0.3,0.0001);
			file<<"Energy\n";
			for(double p:fit)file<<p<<" ";
			file<<"\n";
		}
		{//Energy
			printf("\nEnergy...\n");
			auto data=ReadWeightedFrom2D(
				simulation+"He3.E.FTH1.simulation.txt",
				0.0,0.5,100,
				0.0,0.5,100,
				[](double&,double&){return true;}
			);
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
			Plot<double>().Line("He3.E.FTH1.calibration",[&fit](double x){return fit(ParamSet(x));},0.0,0.3,0.0001);
			file<<"Energy\n";
			for(double p:fit)file<<p<<" ";
			file<<"\n";
		}
		{//theta
			printf("\ntheta...\n");
			auto data=ReadWeightedFrom2D(
				simulation+"He3.th.simulation.txt",
				0.0,0.2,200,
				0.0,0.2,200,
				[](double&x,double&y){return pow(x-y-0.01,2)<0.0003;}
			);
			printf("%i points with weights\n",data->count());
			typedef PolynomFunc<0,0,2> fnc;
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
			Plot<double>().Line("He3.th.calibration",[&fit](double x){return fit(ParamSet(x));},0.0,0.15,0.0001);
			file<<"Theta\n";
			for(double p:fit)file<<p<<" ";
			file<<"\n";
		}
		{//phi
			printf("\nphi...\n");
			auto data=ReadWeightedFrom2D(
				simulation+"He3.phi.simulation.txt",
				0.0,6.3,50,
				0.0,6.3,50,
				[](double&x,double&y){return pow(x-y,2)<0.5;}
			);
			printf("%i points with weights\n",data->count());
			typedef PolynomFunc<0,0,2> fnc;
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
			Plot<double>().Line("He3.phi.calibration",[&fit](double x){return fit(ParamSet(x));},0.0,6.283,0.0005);
			file<<"Phi\n";
			for(double p:fit)file<<p<<" ";
			file<<"\n";
		}
		file.close();
	}
}
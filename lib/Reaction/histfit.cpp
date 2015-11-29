// this file is distributed under 
// GPL v 3.0 license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <functional>
#include <fit.h>
#include <paramfunc.h>
#include <filter.h>
#include <initialconditions.h>
#include <equation.h>
#include <gethist.h>
using namespace std;
using namespace Genetic;
hist FitHistByHists(
	const hist&data,const std::vector<hist>&MC,Genetic::RANDOM&engine,
	const std::vector<double*> out,const std::vector<double*> out_err
){
	auto histsum=[&data,&MC](const ParamSet&P){
		hist sum=data.CloneEmptyBins();
		for(size_t i=0;i<MC.size();i++){
			hist tmp=MC[i];
			tmp*=P[i];
			sum+=tmp;
		}
		return sum;
	};
	printf("Fitting...\n");
	SearchMin<DifferentialMutations<Parabolic>> fit([&data,&MC,histsum](const ParamSet&P){
		hist bg=histsum(P);
		return data.HowClose(bg);
	});
	fit.SetFilter([](const ParamSet&P){
		bool res=true;
		for(double p:P)res&=(p>0);
		return res;
	});
	printf("Fitting init...\n");
	auto init=make_shared<GenerateByGauss>();
	for(auto H:MC)init<<make_pair(1,1);
	fit.Init(MC.size()*15,init,engine);
	printf("Fitting loop...\n");
	while(!fit.RelativeOptimalityExitCondition(0.00000001))
		fit.Iterate(engine);
	printf("Fitting done.\n");
	for(int i=0;(i<fit.Parameters().Count())&&(i<out.size());i++)
		*(out[i])=fit[i];
	for(int i=0;(i<fit.Parameters().Count())&&(i<out_err.size());i++){
		*(out_err[i])=fit.GetParamParabolicError(0.0000001,i);
	}
	return histsum(fit.Parameters());
}
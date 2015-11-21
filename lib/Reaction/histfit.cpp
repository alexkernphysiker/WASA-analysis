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
void AnalyseMMSpectra(hist::point&out_bin,const hist&data,const vector<hist>&MC, RANDOM&engine,function<void(hist&)>handle){
	auto histsum=[&data,&MC](const ParamSet&P){
		hist sum=data.CloneEmptyBins();
		for(size_t i=0;i<MC.size();i++){
			hist tmp=MC[i];
			tmp*=P[i];
			sum+=tmp;
		}
		return sum;
	};
	SearchMin<DifferentialMutations<Parabolic>> fit([&data,&MC,histsum](const ParamSet&P){
		hist bg=histsum(P);
		return data.HowClose(bg);
	});
	fit.SetFilter([](const ParamSet&P){
		bool res=true;
		for(double p:P)res&=(p>0);
		return res;
	});
	fit.SetThreadCount(1);
	auto init=make_shared<GenerateByGauss>();
	for(auto H:MC)init<<make_pair(1,1);
	fit.Init(MC.size()*30,init,engine);
	while(!fit.AbsoluteOptimalityExitCondition(0.0000001))
		fit.Iterate(engine);
	hist fithist=histsum(fit.Parameters());
	handle(fithist);
	out_bin.y=fit[0];
	out_bin.dy=fit.GetParamParabolicError(0.0000001,0);
}
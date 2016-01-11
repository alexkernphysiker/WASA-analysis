// this file is distributed under 
// MIT license
#ifndef JUDIQVAJ
# define JUDIQVAJ
#include <iostream>
#include <string>
#include <sstream>
#include <Genetic/fit.h>
#include <plotstream.h>
namespace SimulationDataProcess{
	using namespace std;
	using namespace Genetic;
	using namespace PlotStream;
	string SimulationDataPath();
	shared_ptr<FitPoints> WeightPoints(shared_ptr<FitPoints>src,ParamSet&&pos={},ParamSet&&delta={},size_t index=0);
	template<class FITFUNC>
	void ProcessFit(
		string&&reconstructionname,
		shared_ptr<IInitialConditions>init,
		shared_ptr<IParamCheck>filter,
		RANDOM&R
	){
		auto points=make_shared<FitPoints>();
		ifstream file;
		file.open(SimulationDataPath()+reconstructionname+".simulation.txt");
		if(file){
			cout<<"reading...";
			string line;
			while(getline(file,line)){
				istringstream str(line);
				ParamSet X;
				str>>X;
				double y;
				X>>y;
				points<<Point(X,y,1);
			}
			file.close();
			cout<<"done."<<endl;
		}
		cout<<"Init1"<<endl;
		auto weighted=WeightPoints(points);
		Binner weightedpic(
			reconstructionname+"_weighted",
			Binner::binparam(1,0,0.3,10)
		);
		weightedpic.Fill([](size_t i,string&&name){
			return make_shared<SimplePlotStream>(name+"_"+to_string(i),make_pair(0,2));
		});
		for(Point&P:*weighted)
			if(P.wy()>=1)
				weightedpic<<P;
		/*FitFunction<DifferentialMutations<>,FITFUNC,SumWeightedSquareDiff> fit(weighted);
		cout<<"Init2"<<endl;
		fit.SetFilter(filter).Init(30*FITFUNC::ParamCount,init,R);
		cout<<"Fitting"<<endl;
		while(!fit.AbsoluteOptimalityExitCondition(0.000001)){
			fit.Iterate(R);
			cout<<fit.Optimality()<<"<S<"<<fit.Optimality(fit.PopulationSize()-1)<<"     \r";
		}
		cout<<endl;
		cout<<"Fit parameters:"<<endl<<fit.Parameters()<<endl;
		cout<<"Errors:"<<endl<<fit.GetParamParabolicErrors(parEq(fit.ParamCount(),0.001))<<endl;
		ofstream out;
		out.open(reconstructionname+".fit.txt");
		if(out){
			out<<fit;
			out.close();
		}*/
	}
};
#endif
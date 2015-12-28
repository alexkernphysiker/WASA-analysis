// this file is distributed under 
// MIT license
#ifndef JUDIQVAJ
# define JUDIQVAJ
#include <iostream>
#include <string>
#include <sstream>
#include <Genetic/fit.h>
std::string SimulationDataPath();
std::shared_ptr<Genetic::FitPoints> WeightPoints(
	std::shared_ptr<Genetic::FitPoints>src,
	Genetic::ParamSet&&pos={},Genetic::ParamSet&&delta={},
	size_t index=0
);
template<class FITFUNC>
void ProcessFit(
	std::string reconstructionname,
	std::shared_ptr<Genetic::IInitialConditions>init,
	std::shared_ptr<Genetic::IParamCheck>filter,
	Genetic::RANDOM&R
){
	using namespace std;
	using namespace Genetic;
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
	FitFunction<DifferentialMutations<>,FITFUNC,SumWeightedSquareDiff> fit(WeightPoints(points));
	cout<<"Init2"<<endl;
	fit.SetFilter(filter).Init(20*FITFUNC::ParamCount,init,R);
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
	}
}
#endif
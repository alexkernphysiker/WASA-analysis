// this file is distributed under 
// MIT license
#ifndef JUDIQVAJ
# define JUDIQVAJ
#include <iostream>
#include <string>
#include <sstream>
#include <Genetic/fit.h>
#include <Genetic/paramsort.h>
namespace SimulationDataProcess{
	using namespace std;
	using namespace Genetic;
	string SimulationDataPath();
	template<class FITFUNC>
	void ProcessEnergyThetaFit(
		string&&reconstructionname,pair<double,double> E_range,
		shared_ptr<IInitialConditions>init,
		shared_ptr<IParamCheck>filter,
		RANDOM&R
	){
		auto params_shown=make_pair(0,2);
		auto theta_binning=BinningParam(1,make_pair(0.1,0.16),10);
		auto Edep_binning=BinningParam(0,E_range,100);
		auto Ek_binning=BinningParam(2,E_range,100);
		ParamsPerBinsCounter<3> Binner({theta_binning,Edep_binning,Ek_binning});
		auto AllData=make_shared<FitPoints>();
		ifstream file;
		file.open(SimulationDataPath()+reconstructionname+".simulation.txt");
		if(file){
			cout<<"reading..."<<endl;
			string line;
			while(getline(file,line)){
				istringstream str(line);
				ParamSet X;
				str>>X;
				double y;
				X>>y;
				Binner<<X;
				AllData<<Point(X,y);
			}
			file.close();
			cout<<"done."<<endl;
		}
		cout<<"Init1"<<endl;
		auto points=make_shared<FitPoints>();
		Binner.FullCycle([&points](ParamSet&P,unsigned long cnt){
			double y;
			ParamSet X=P;X>>y;
			points<<Point(X,y,cnt);
		});
		FitFunction<DifferentialMutations<>,FITFUNC,SumWeightedSquareDiff> fit(points);
		cout<<"Init2"<<endl;
		fit.SetFilter(filter).Init(30*FITFUNC::ParamCount,init,R);
		cout<<"Fitting"<<endl;
		while(
			(!fit.AbsoluteOptimalityExitCondition(0.000001))&&
			(!fit.RelativeOptimalityExitCondition(0.001))
		){
			fit.Iterate(R);
			cout<<fit.Optimality()<<"<S<"<<fit.Optimality(fit.PopulationSize()-1)<<"     \r";
		}
		cout<<endl;
		cout<<"Fit parameters:"<<endl<<fit.Parameters()<<endl;
		{
			ofstream out;
			out.open(reconstructionname+".fit.txt");
			if(out){
				out<<fit;
				out.close();
			}
		}
		SimplePlotStream total_pic(reconstructionname+"_total",params_shown);
		for(Point&P:*AllData){
			total_pic<<(ParamSet(P.X())<<P.y());
		}
	}
};
#endif
// this file is distributed under 
// MIT license
#ifndef JUDIQVAJ
# define JUDIQVAJ
#include <iostream>
#include <string>
#include <sstream>
#include <math_h/error.h>
#include <Genetic/fit.h>
#include <Genetic/paramsort.h>
namespace SimulationDataProcess{
	using namespace std;
	using namespace Genetic;
	using namespace MathTemplates;
	string SimulationDataPath();
	template<class FITFUNC>
	void ProcessEnergyThetaFit(
		string&&reconstructionname,pair<double,double> E_range,
		shared_ptr<IInitialConditions>init,
		shared_ptr<IParamCheck>filter,
		RANDOM&R
	){
		auto params_shown=make_pair(0,2);
		auto theta_binning=BinningParam(1,make_pair(0.1,0.16),12);
		auto Edep_binning=BinningParam(0,E_range,25);
		auto Ek_binning=BinningParam(2,E_range,25);
		
		ParamsPerBins forplots(theta_binning);
		ParamsPerBinsCounter<3> Binner({theta_binning,Edep_binning,Ek_binning});
		auto AllData=make_shared<FitPoints>();{
			ifstream file;
			file.open(SimulationDataPath()+reconstructionname+".simulation.txt");
			if(file){
				cout<<"reading..."<<endl;
				string line;
				while(getline(file,line)){
					istringstream str(line);
					ParamSet X;
					str>>X;
					Binner<<X;
					forplots<<X;
					double y;
					X>>y;
					AllData<<Point(X,y);
				}
				file.close();
				cout<<"done."<<endl;
			}else
				throw Exception<ifstream>("No input data");
		}
		cout<<"Init1"<<endl;
		auto points=make_shared<FitPoints>();
		Binner.FullCycle([&points](const ParamSet&center_pos,const unsigned long events_count){
			double y;
			ParamSet X=center_pos;X>>y;
			points<<Point(X,y,events_count);
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
		cout<<"done.                                                                            "<<endl;
		{
			ofstream out;
			out.open(reconstructionname+".fit.txt");
			if(out){
				out<<fit;
				out.close();
			}else
				throw Exception<ofstream>("Cannot write output");
		}
		SimplePlotStream total_pic(reconstructionname+" total",params_shown);
		for(Point&P:*AllData){
			total_pic<<(ParamSet(P.X())<<P.y());
		}
		for(size_t i=0;i<forplots.count();i++){
			double theta=forplots.bin_center(i);
			SimplePlotStream plot(reconstructionname+" theta="+to_string(theta),params_shown);
			for(ParamSet&P:forplots[i])plot<<P;
			plot.AddFunc([theta,&fit](double Edep){return fit({Edep,theta});});
		}
	}
};
#endif
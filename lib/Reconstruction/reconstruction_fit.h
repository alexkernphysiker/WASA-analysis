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
		auto Edep_binning=BinningParam(0,E_range,30);
		auto theta_binning=BinningParam(1,make_pair(0.10,0.16),12);
		auto Ek_binning=BinningParam(2,E_range,30);
		
		ParamsPerBinsCounter<3> Binner({Edep_binning,theta_binning,Ek_binning});
		auto AllData=vector<ParamSet>();
		{
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
					AllData.push_back(X);
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
		fit.SetFilter(filter).Init(25*FITFUNC::ParamCount,init,R);
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
		ParamsPerBins binned_plots(theta_binning);
		{
			SimplePlotStream total_pic(reconstructionname+" weights",params_shown);
			for(const Point&P:*points)if(P.wy()>5){
				total_pic<<(ParamSet(P.X())<<P.y());
				//binned_plots<<(ParamSet(P.X())<<P.y());
			}
			SimplePlotStream total_points(reconstructionname+" points",params_shown);
			for(const ParamSet&P:AllData){
				total_points<<P;
				binned_plots<<P;
			}
		}
		for(size_t i=0;i<binned_plots.count();i++){
			double theta=binned_plots.bin_center(i);
			SimplePlotStream plot(reconstructionname+" theta="+to_string(theta),params_shown);
			for(ParamSet&P:binned_plots[i])plot<<P;
			plot.AddFunc([theta,&fit](double Edep){return fit({Edep,theta});});
		}
	}
};
#endif
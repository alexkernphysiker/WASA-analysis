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
	template<class FITFUNC>
	void ProcessEnergyThetaFit(
		string&&reconstructionname,pair<double,double> E_range,
		shared_ptr<IInitialConditions>init,
		shared_ptr<IParamCheck>filter,
		RANDOM&R
	){
		auto params_shown=make_pair(0,2);
		auto theta_binning=Binner::binparam(1,0.1,0.16,10);
		auto substream=[&params_shown](size_t i,string&&name){
			return make_shared<SimplePlotStream>(name+"_"+to_string(i),params_shown);
		};
		auto AllData=make_shared<FitPoints>();
		auto points=make_shared<FitPoints>();
		double de=0.005,dth=0.005,dek=0.005;
		for(double e=E_range.first;e<E_range.second;e+=(de*2.0))
			for(double th=0.1;th<0.16;th+=(dth*2.0))
				for(double ek=E_range.first;ek<E_range.second;ek+=(dek*2.0))
					points<<Point({e,th},ek,0);
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
				AllData<<Point(X,y);
				for(Point&P:*points){
					if((abs(P.y()-y)<dek)&&(abs(P.X()[0]-X[0])<de)&&(abs(P.X()[1]-X[1])<dth))
						P.WY()+=1;
				}  
			}
			file.close();
			cout<<"done."<<endl;
		}
		cout<<"Init1"<<endl;
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
		cout<<"Errors:"<<endl<<fit.GetParamParabolicErrors(parEq(fit.ParamCount(),0.001))<<endl;
		{
			ofstream out;
			out.open(reconstructionname+".fit.txt");
			if(out){
				out<<fit;
				out.close();
			}
		}
		auto plotfunc=[&fit](double theta,double E){return fit({E,theta});};
		Binner bined_pic(reconstructionname+"_bined",theta_binning);
		bined_pic.Fill(substream).AddFunc(plotfunc);
		for(Point&P:*points)if(P.wy()>2)bined_pic<<P;
		
		SimplePlotStream total_pic(reconstructionname+"_total",params_shown);
		Binner bined_data(reconstructionname+"_data",theta_binning);
		bined_data.Fill(substream).AddFunc(plotfunc);
		for(Point&P:*AllData){
			bined_data<<P;
			total_pic<<P;
		}
	}
};
#endif
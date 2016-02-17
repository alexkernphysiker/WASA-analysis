// this file is distributed under 
// MIT license
#ifndef JUDIQVAJ
# define JUDIQVAJ
#include <iostream>
#include <string>
#include <sstream>
#include <gnuplot_wrap.h>
#include <math_h/error.h>
#include <math_h/hist.h>
#include <Genetic/fit.h>
namespace SimulationDataProcess{
	using namespace std;
	using namespace Genetic;
	using namespace MathTemplates;
	using namespace GnuplotWrap;
	string SimulationDataPath();
	template<class FITFUNC>
	void He3ForEtaFit(string&&reconstructionname,pair<double,double>&&E_range,shared_ptr<IInitialConditions>init,RANDOM&R){
		auto params_shown=make_pair(0,2);
		vector<value<double>> theta_bins=BinsByStep(0.10,0.002,0.13);
		vector<Distribution2D<double>> E_sp2;
		for(const value<double>&bin:theta_bins)
			E_sp2.push_back(Distribution2D<double>(BinsByStep(E_range.first,0.005,E_range.second),BinsByStep(E_range.first,0.005,E_range.second)));
		cout<<theta_bins.size()<<" theta bins"<<endl;
		{
			ifstream file;
			file.open(SimulationDataPath()+reconstructionname+".simulation.txt");
			if(file){
				cout<<"reading..."<<endl;
				string line;
				while(getline(file,line)){
					istringstream str(line);
					double Edep,theta,Ekin;
					str>>Edep>>theta>>Ekin;
					for(size_t i=0,n=theta_bins.size();i<n;i++)
						if(theta_bins[i].contains(theta))
							E_sp2[i]<<make_pair(Edep,Ekin);
				}
				file.close();
				cout<<"done."<<endl;
			}else
				throw Exception<ifstream>("No input data");
		}
		cout<<"Init1"<<endl;
		auto points=make_shared<FitPoints>();
		for(size_t i=0,n=theta_bins.size();i<n;i++)
			E_sp2[i].FullCycle([&theta_bins,&E_sp2,&points,i](Distribution2D<double>::Point&&P){
				if(!P.Z().contains(0))
					points<<Point({P.X().val(),theta_bins[i].val()},P.Y().val(),P.Z().val());
			});
		cout<<points->size()<<" points"<<endl;
		FitFunction<DifferentialMutations<>,FITFUNC,SumWeightedSquareDiff> fit(points);
		cout<<"Init2"<<endl;
		fit.Init(25*FITFUNC::ParamCount,init,R);
		cout<<"population "<<fit.PopulationSize()<<endl;
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
		for(size_t i=0,n=theta_bins.size();i<n;i++){
			PlotDistribution2D<double>(sp2).Distr(E_sp2[i],"Theta=["+to_string(theta_bins[i].min())+":"+to_string(theta_bins[i].max())+"]");
			typedef PlotDistribution2D<double>::point P;
			vector<P> func_line;
			for(const value<double>Ed:BinsByStep(E_range.first,0.001,E_range.second))
				func_line.push_back(P(Ed.val(),fit({Ed.val(),theta_bins[i].val()}),500));
			PlotDistribution2D<double>(normal).Distr(E_sp2[i],"Theta=["+to_string(theta_bins[i].min())+":"+to_string(theta_bins[i].max())+"]")
				.Line(func_line,"Theta="+to_string(theta_bins[i].val()));
		}
	}
};
#endif
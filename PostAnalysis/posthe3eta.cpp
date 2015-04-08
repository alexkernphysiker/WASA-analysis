#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <unistd.h>

#include <paramfunc.h>
#include <fitpoints.h>
#include <filter.h>
#include <initialconditions.h>
#include <genetic.h>

#include <phys_constants.h>
#include "gethist.h"
#include "gnuplot.h"
using namespace std;
using namespace Fit;
typedef Mul<Par<0>,Func3<Gaussian,Arg<0>,Par<1>,Par<2>>> Foreground;
int main(int,char**){
#include "env.cc"
	string MCFile=inputpath+"/MCHe3Eta.root";
	outpath+="/He3Eta";
	SetPlotOutput(outpath);
	vector<string> Kinematics;
	Kinematics.push_back("Histograms");
	Kinematics.push_back("Kinematics");
	vector<string> Reconstruction;
	Reconstruction.push_back("Histograms");
	Reconstruction.push_back("Reconstruction");
	vector<string> EventsCount;
	EventsCount.push_back("Histograms");
	EventsCount.push_back("OutputEventsCount");
	LinearInterpolation<double> mc_norm,mc_dnorm;{
		hist normhist(MCFile,Reconstruction,"P_beam");
		for(hist::point p:normhist){
			mc_norm<<make_pair(p.x,p.y);
			mc_dnorm<<make_pair(p.x,p.dy);
		}
	}
	vector<pair<double,shared_ptr<FitPoints>>> missingmass;{
		LinearInterpolation<double> acceptance,dacceptance;{
			hist registered(MCFile,EventsCount,"DependenceOnBeam");
			for(hist::point histpoint:registered){
				double norm=mc_norm(histpoint.x);
				double dnorm=mc_dnorm(histpoint.x);
				if(mc_norm(histpoint.x)>50000){
					string subhistname=string("MissingMass")+to_string(int(histpoint.x*1000));
					hist subhist(MCFile,Kinematics,subhistname);
					auto points=make_shared<FitPoints>();
					for(hist::point subhistpoint:subhist)
						if((subhistpoint.x>=0.53)&&(subhistpoint.x<=0.56)){
							FitPoints::Point point;
							point.X<<subhistpoint.x;
							point.WX<<subhistpoint.dx;
							point.y=subhistpoint.y/subhistpoint.dx;
							point.wy=((subhistpoint.dy>1.0)?subhistpoint.dy:1.0)/subhistpoint.dx;
							points<<point;
						}
					missingmass.push_back(make_pair(histpoint.x,points));
					acceptance<<make_pair(histpoint.x,histpoint.y/norm);
					dacceptance<<make_pair(histpoint.x,(dnorm*histpoint.y/pow(norm,2))+(histpoint.dy/norm));
				}
			}
		}
		Plot fig1;
		fig1.Points("Estimated Acceptance",acceptance,[&dacceptance](double x){return dacceptance(x);});
	}
	LinearInterpolation<double> mc_acc,mc_dacc;
	for(auto hist:missingmass){
		double pbeam=hist.first;
		string ind=to_string(int(pbeam*1000));
		Plot fit_image;
		fit_image.Points(ind+"points",hist.second);
		printf("Fitting missing mass histogram for p=%f GeV/c...\n",hist.first);
		DifferentialRandomMutations<> fit(make_shared<Foreground>(),ChiSquareWithXError(hist.second));
		fit.SetFilter(make_shared<And>()
			<<(make_shared<Above>()<<0<<0.5<<0)
			<<(make_shared<Below>()<<INFINITY<<0.6<<0.2)
			<<make_shared<Filter<>>([](ParamSet &P){return pow(P[1]-m_eta,2)<0.5;})
		);
		double norm=mc_norm(hist.first);
		double dnorm=mc_dnorm(hist.first);
		printf("norm=%f\n",norm);
		fit.Init(30,make_shared<GenerateByGauss>()
			<<make_pair(2.0*norm/3.0,norm/2.0)
			<<make_pair(m_eta,0.01)
			<<make_pair(0.05,0.05)
		);
		printf("Inited...\n");
		while(!fit.AbsoluteOptimalityExitCondition(0.0000001)){
			fit.Iterate();
			printf("%i iterations; %f<=chi^2<=%f         \r",
				fit.iteration_count(),
				fit.Optimality(),
				fit.Optimality(fit.PopulationSize()-1)
			);
		}
		printf("\nParameters:\n");
		for(double p:fit)
			printf("\t%f;",p);
		printf("\nErrors:\n");
		ParamSet err=fit.GetParamParabolicError(ParamSet(1,0.001,0.001));
		for(double p:err)
			printf("\t%f;",p);
		mc_acc<<make_pair(pbeam,fit[0]/norm);
		mc_dacc<<make_pair(pbeam,(dnorm*fit[0]/pow(norm,2))+(err[0]/norm));
		printf("\n>>>>> %f+/-%f\n",mc_acc(pbeam),mc_dacc(pbeam));
		fit_image.Function(ind+"fit",[&fit](double x){return fit(ParamSet(x));},0.5,0.6,0.0001);
	}
	Plot fig_acc;
	fig_acc.Points("Acceptance",mc_acc,[&mc_dacc](double x){return mc_dacc(x);});
}
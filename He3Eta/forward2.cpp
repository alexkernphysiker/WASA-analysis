// this file is distributed under 
// GPL license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Genetic/fit.h>
#include <Genetic/paramfunc.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Genetic/parabolic.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include "he3eta.h"
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main(){
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_luminosity");
    vector<string> histpath_forward_reconstr={"Histograms","He3Forward_Reconstruction"};
    hist<double> norm=Hist(MC,"He3eta",histpath_forward_reconstr,"0-Reference");
    hist<double> luminosity,luminosity2,data_chi_sq;
    vector<hist<double>> parhists;
    RANDOM r_eng;
    for(size_t bin_num=0,bin_count=norm.size();bin_num<bin_count;bin_num++)
	if(norm[bin_num].X()>5.0){
	    const auto&Q=norm[bin_num].X();
	    const auto&N=norm[bin_num].Y();
	    string Qmsg="Q in ["+to_string(Q.min())+":"+to_string(Q.max())+"] MeV";
	    const hist<double> data=Hist(DATA,"",histpath_forward_reconstr,
		string("MissingMass-Bin-")+to_string(bin_num)
	    ).XRange(0.525,0.57);
	    const auto chain=ChainWithStep(0.525,0.001,0.57);
	    const hist<double> mc=Hist(MC,"He3eta",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num))/N;
	    const LinearInterpolation<double> fg=mc.toLine();
	    const auto&data_count=data.TotalSum().val();
	    auto BG=[&data_count](const ParamSet&X,const ParamSet&P){
		const double res=data_count*Polynom(X[0],P,3,2);
		return (res>0)?res:0.0;
	    };
	    Fit<AbsoluteMutations<DifferentialMutations<Uncertainty>>>
	    FIT(make_shared<FitPoints>(data),
		[&fg,BG](const ParamSet&X,const ParamSet&P){
		    return P[0]*fg(X[0]+P[1])+BG(X,P);
		}
	    );
	    FIT.SetFilter([BG](const ParamSet&P){
		return (P[0]>0)&&(BG({-P[3]/(2.0*P[4])},P)>0)
		    &&(P[2]<0)&&(P[3]>0)&&(P[4]<0)&&(P[5]<0);
	    });
	    auto init=make_shared<InitialDistributions>()
		<<make_shared<DistribUniform>(0.0,data_count*Q.val()/100.)
		<<make_shared<DistribGauss>(0.,0.0001)
		<<make_shared<DistribGauss>(-10,10.)
		<<make_shared<DistribGauss>( 10,10)
		<<make_shared<DistribGauss>(-5,5)
		<<make_shared<DistribGauss>(-3,2)
	    ;
	    FIT.SetAbsoluteMutationCoefficients({1.,0.0001,1.,1.,1.,1.});
	    FIT.SetAbsoluteMutationsProbability(0.2);
	    FIT.SetUncertaintyCalcDeltas({0.1,0.0001,0.01,0.01,0.01,0.001});
	    FIT.Init(500,init,r_eng);

	    cout<<endl;
	    while(
		!FIT.AbsoluteOptimalityExitCondition(0.0000001)
	    ){
		FIT.Iterate(r_eng);
		cout<<"Fitting: "<<FIT.iteration_count()<<" iterations; "
		<<FIT.Optimality()<<"<chi^2<"
		<<FIT.Optimality(FIT.PopulationSize()-1)
		<<"          \r";
	    }
	    const auto&P=FIT.ParametersWithUncertainties();
	    if(parhists.size()==0){
		for(size_t i=0;i<P.size();i++)
		    parhists.push_back(hist<double>());
	    }
	    for(size_t i=0;i<P.size();i++)
		parhists[i]<< point<value<double>>(Q,P[i]);
	    data_chi_sq << point<value<double>>(Q,FIT.Optimality()/(data.size()-FIT.ParamCount()));
	    cout<<endl;
	    Plot<double> exp_plot;
	    exp_plot.Hist(data,"DATA")
	    << "set key on"<< "set title '"+Qmsg+"'"
	    << "set xlabel 'Missing mass, GeV'"
	    << "set ylabel 'counts'"
	    << "set yrange [-200:]"<<"unset log y";
	    const SortedPoints<double>
		totalfit([&FIT](double x)->double{return FIT({x});},chain),
		background([&FIT,BG](double x)->double{return BG({x},FIT.Parameters());},chain);
	    exp_plot.Line(totalfit,"fit").Line(background,"background");

	    hist<double> bg;
	    for(const auto&po:data){
		const double&x=po.X().val();
		const ParamSet& p=FIT.Parameters();
		double v=BG({x},p),u=0;
		for(size_t index=1;index<p.size();index++){
		    ParamSet p1=p,p2=p1;
		    const double&delta=P[index].uncertainty();
		    p1(index)+=delta;
		    p2(index)-=delta;
		    u+=pow(BG({x},p1)-BG({x},p2),2);
		}
		bg<<point<value<double>>(po.X(),{v,sqrt(u)});
	    }
	    hist<double> clean=data-bg;
	    Plot<double> subplot;
	    subplot.Object("0 title \"\"").Hist(clean);
	    subplot.Hist(clean=clean.XRange(0.540,0.556))
	    << "set key on"<< "set title '"+Qmsg+"'"
	    << "set xlabel 'Missing mass, GeV'"
	    << "set ylabel 'counts'"
	    << "set yrange [-200:]"<<"unset log y";

	    luminosity << point<value<double>>(Q,
		(P[0]/he3eta_sigma()(Q))
		*double(trigger_he3_forward.scaling)
	    );
	    luminosity2 << point<value<double>>(Q,
		((clean.TotalSum()/mc.TotalSum())/he3eta_sigma()(Q))
		*double(trigger_he3_forward.scaling)
	    );
	}
    for(size_t i=0;i<parhists.size();i++)
	Plot<double>().Hist(parhists[i])
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'parameter"+to_string(i)+"'";

    Plot<double>().Hist(data_chi_sq)
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'chi^2/d, n.d.'" 
    << "set yrange [0:]"<<"unset log y";

    Plot<double>()
    .Hist(hist<double>(he3eta_sigma().func(),BinsByStep(5.0,2.5,30.0)))
    .Hist(he3eta_sigma(),"Data from other experiments")
    << "set title 'Cross section of He3eta used in the calculations'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'sigma(^3He eta), nb'"
    << "set xrange [0:45]"<< "set yrange [0:600]";

    auto runs=PresentRuns("");
    Plot<double>().Hist(luminosity2,"Substracted").Hist(luminosity,"Fit parameter")
    << "set title 'Integral luminosity estimation ("+to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs)'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Integral luminosity, nb^{-1}'" 
    << "set xrange [0:45]"<< "set yrange [0:]";
    Plot<double>().Hist(luminosity2)
    << "set title 'Integral luminosity estimation ("+to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs)'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Integral luminosity, nb^{-1}'" 
    << "set xrange [0:45]"<< "set yrange [0:]";
}


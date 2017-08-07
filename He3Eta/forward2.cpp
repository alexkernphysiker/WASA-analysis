// this file is distributed under 
// GPL license
#include <iostream>
#include <iomanip>
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
    const auto runs=PresentRuns("F");
    const string runmsg=to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs";
    hist<double> norm=Hist(MC,"He3eta",histpath_forward_reconstr,"0-Reference");
    hist<double> luminosity,data_chi_sq;
    vector<hist<double>> parhists;
    RANDOM r_eng;
    for(size_t bin_num=0,bin_count=norm.size();bin_num<bin_count;bin_num++)
	if(norm[bin_num].X()>2.5){
	    const auto&Q=norm[bin_num].X();
	    const auto&N=norm[bin_num].Y();
	    const string Qmsg=static_cast<stringstream&>(stringstream()
		<<"Q in ["<<setprecision(3)
		<<Q.min()<<"; "<<Q.max()<<"] MeV"
	    ).str();
	    const hist<double> data=Hist(DATA,"F",histpath_forward_reconstr,
		string("MissingMass-Bin-")+to_string(bin_num)
	    ).XRange(0.53,0.57);
	    const auto chain=ChainWithStep(0.53,0.0001,0.57);
	    const auto cut=make_pair(0.541,0.554);
	    const hist<double> mc=Hist(MC,"He3eta",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num))/N;
	    const LinearInterpolation<double> fg=mc.toLine();
	    const auto&data_count=data.TotalSum().val();
	    const auto BG=[&data_count](const ParamSet&X,const ParamSet&P){
		const double res=data_count*Polynom(X[0],P,3,0);
		return (res>0)?res:0.0;
	    };
	    const auto peak_reg= data.XRange(cut.first,cut.second);
	    const auto data_bg= data.XExclude(cut.first,cut.second);
	    cout<<endl<<Qmsg<<endl<<endl;
	    Fit<AbsoluteMutations<DifferentialMutations<Uncertainty>>> FIT(make_shared<FitPoints>(data_bg),BG);
	    FIT
	    .SetAbsoluteMutationCoefficients({1.0,1.0,1.0,1.0})
	    .SetAbsoluteMutationsProbability(0.2)
	    .SetUncertaintyCalcDeltas({0.1,0.1,0.1,0.1})
	    .SetFilter([BG,&peak_reg](const ParamSet&P){
		bool valid=(BG({peak_reg.left().X().val()},P)>0.0);
		if(valid)for(const auto&p:peak_reg){
		    valid&=(p.Y().max()>BG({p.X().val()},P));
		}
		return valid;
	    });
	    auto init=make_shared<InitialDistributions>()
		<<make_shared<DistribGauss>(0,5)
		<<make_shared<DistribGauss>(0,5)
		<<make_shared<DistribGauss>(0,5)
		<<make_shared<DistribGauss>(0,5)
	    ;
	    FIT.Init(300,init,r_eng);
	    cout<<endl;
	    while(!FIT.AbsoluteOptimalityExitCondition(0.0000001)){
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
	    exp_plot.Hist(data).Hist(data_bg)
	    << "set key on"<< "set title '"+Qmsg+", "+runmsg+"'"
	    << "set xlabel 'Missing mass, GeV'"
	    << "set ylabel 'counts'"
	    << "set yrange [-200:]"<<"unset log y";
	    const SortedPoints<double>
		background([&FIT,BG](double x)->double{return BG({x},FIT.Parameters());},chain);
	    exp_plot.Line(background,"background");
	    hist<double> bg;
	    for(const auto&po:data){
		const double&x=po.X().val();
		const ParamSet& p=FIT.Parameters();
		double v=BG({x},p),u=0;
		for(size_t index=0;index<p.size();index++){
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
	    subplot.Hist(clean);
	    subplot.Hist(clean=clean.XRange(cut.first,cut.second)).Object("0 title \"\"")
	    << "set key on"<< "set title '"+Qmsg+", "+runmsg+"'"
	    << "set xlabel 'Missing mass, GeV'"
	    << "set ylabel 'counts'"
	    << "set yrange [-200:]"<<"unset log y";

	    luminosity << point<value<double>>(Q,
		((clean.TotalSum()/mc.TotalSum())/he3eta_sigma()(Q))
		* double(trigger_he3_forward.scaling)
	    );
	}
    for(size_t i=0;i<parhists.size();i++)
	Plot<double>().Hist(parhists[i])
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'parameter"+to_string(i)+"'";
    Plot<double>().Hist(data_chi_sq)
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'chi^2/d, n.d.'" 
    << "set yrange [0:3]"<<"unset log y";

    Plot<double>()
    .Hist(hist<double>(he3eta_sigma().func(),BinsByStep(2.5,2.5,30.0)))
    .Hist(he3eta_sigma(),"Data from other experiments")
    << "set title 'Cross section of He3eta used in the calculations'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'sigma(^3He eta), nb'"
    << "set xrange [0:45]"<< "set yrange [0:600]";

    Plot<double>().Hist(luminosity)
    << "set title 'Integrated luminosity estimation, "+runmsg+"'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Integrated luminosity, nb^{-1}'" 
    << "set xrange [0:30]"<< "set yrange [0:]";

    Plot<double>().Hist(luminosity*runs.second/runs.first)
    << "set title 'Integrated luminosity estimation, all runs'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Integrated luminosity, nb^{-1}'" 
    << "set xrange [0:30]"<< "set yrange [0:]";
}


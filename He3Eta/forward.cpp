// this file is distributed under 
// GPL license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Genetic/searchmin.h>
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
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"he3eta_forward");
    vector<string> histpath_forward_reconstr={"Histograms","He3Forward_Reconstruction"};
    vector<string> reaction={"He3eta","He3pi0pi0pi0","He3pi0pi0","He3pi0"};
    vector<hist<double>> norm;
    for(const string& r:reaction)
	norm.push_back(Hist(MC,r,histpath_forward_reconstr,"0-Reference"));

    vector<hist<double>> acceptance;
    for(size_t i=0;i<norm.size();i++)
	acceptance.push_back(hist<double>());
    hist<double> luminosity,bg_chi_sq;
    vector<hist<double>> fit_params;
    for(size_t i=0;i<reaction.size();i++)
	fit_params.push_back(hist<double>());
    RANDOM r_eng;
    for(size_t bin_num=0,bin_count=norm[0].size();bin_num<bin_count;bin_num++)
	if(norm[0][bin_num].X()>2.5){
	    const auto&Q=norm[0][bin_num].X();
	    string Qmsg="Q in ["+to_string(Q.min())+":"+to_string(Q.max())+"] MeV";
	    auto transform=[](hist<double>&h){h=h.XRange(0.42,0.58);};

	    hist<double> data=Hist(DATA,"",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
	    transform(data);
	    Plot<double> exp_plot;
	    exp_plot.Hist(data,"DATA")
	    << "set key on"<< "set title '"+Qmsg+"'"
	    << "set xlabel 'Missing mass, GeV'"
	    << "set ylabel 'counts'"
	    << "set yrange [0:1500]"<<"unset log y";
	    vector<hist<double>> theory;{
		Plot<double> th_plot;
		th_plot<< "set key on"<< "set xlabel 'Missing mass, MeV'"
		<< "set title '"+Qmsg+"'"
		<< "set ylabel 'acceptance density, GeV^{-1}'";
		for(size_t i=0;i<reaction.size();i++){
		    hist<double> react_sim=Hist(MC,reaction[i],histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
		    transform(react_sim);
		    auto N=norm[i][bin_num].Y();
		    th_plot.Line(react_sim.toLine()/(N.val()*react_sim[0].X().uncertainty()*2.),reaction[i]);
		    auto MN=value<double>::std_error(react_sim.TotalSum().val());
		    acceptance[i] << point<value<double>>(Q,MN/N);
		    theory.push_back(react_sim/N);
		}
	    }

	    SearchMin<DifferentialMutations<ParabolicErrorEstimationFromChisq>>
	    fit([&theory,&data](const ParamSet&P){
		double res=0;
		for(size_t i=1,n=data.size();i<n;i++){
		    value<double> exp_p=data[i].Y(),the_p=0;
		    for(size_t j=0,n=theory.size()-1;j<n;j++)
			the_p+=theory[j][i].Y()*P[j];
		    res+=exp_p.NumCompare(the_p);
		}
		return res;
	    });
	    fit.SetMutationCoefficient(0.8);
	    fit.SetFilter([](const ParamSet&P)->bool{return (P[0]>0)&&(P[1]>0)&&(P[2]>0)&&(P[3]>=0);});
	    const auto&data_count=data.TotalSum().val();
	    fit.Init(100,
		make_shared<InitialDistributions>()
		    <<make_shared<DistribUniform>(0.0,2.0*data_count)
		    <<make_shared<DistribUniform>(0.0,2.0*data_count)
		    <<make_shared<DistribUniform>(0.0,2.0*data_count)
		    <<make_shared<FixParam>(0.0)
		,r_eng
	    );
	    while(!fit.AbsoluteOptimalityExitCondition(0.000000000001)){
		fit.Iterate(r_eng);
		cout<<fit.iteration_count()<<" iterations; "
		<<fit.Optimality()<<"<chi^2<"
		<<fit.Optimality(fit.PopulationSize()-1)
		<<"          \r";
	    }
	    fit.SetUncertaintyCalcDeltas({0.1,1.0,1.0,1.0});
	    const auto&P=fit.ParametersWithUncertainties();
	    for(size_t i=0;i<reaction.size();i++)
		fit_params[i]<< point<value<double>>(Q,P[i]);
	    bg_chi_sq << point<value<double>>(Q,fit.Optimality()/(data.size()-fit.ParamCount()));
	    exp_plot
		.Line(hist<double>(theory[0]*P[0]+theory[1]*P[1]+theory[2]*P[2]).toLine(),"Total")
		.Line(hist<double>(theory[1]*P[1]+theory[2]*P[2]+theory[3]*P[3]).toLine(),"Background")
		.Line(hist<double>(theory[1]*P[1]).toLine(),reaction[1])
		.Line(hist<double>(theory[2]*P[2]).toLine(),reaction[2])
		//.Line(hist<double>(theory[3]*P[3]).toLine(),reaction[3])
	    ;
	    luminosity << point<value<double>>(Q,
	       (P[0]/he3eta_sigma()(Q))
	       *double(trigger_he3_forward.scaling)
	    );
	}
    Plot<double>().Hist(bg_chi_sq) 
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'chi^2/d, n.d.'" 
    << "set yrange [0:]"<<"unset log y";

    {//Plot acceptance
	Plot<double> acc,par;
	acc << "set key on" 
	<< "set title 'Acceptance'"
	<< "set yrange [0:1.0]"<<"unset log y"
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'Acceptance, n.d.'";
	par << "set key on" 
	<< "set title 'Fitted coefficients'"
	<< "set yrange [0:]"
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'coefficient, n.d.'";
	for(size_t i=0;i<reaction.size();i++){
	    acc.Hist(acceptance[i],reaction[i]);
	    par.Hist(fit_params[i],reaction[i]);
	}
    }

    Plot<double>()
    .Hist(hist<double>(he3eta_sigma().func(),BinsByStep(10.0,2.5,30.0)))
    .Hist(he3eta_sigma(),"Data from other experiments")
    << "set title 'Cross section of "+reaction[0]+" used in the calculations'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'sigma(^3He eta), nb'"<<"unset log y"
    << "set xrange [0:45]"<< "set yrange [0:600]";

    const auto runs=PresentRuns("");
    Plot<double>().Hist(luminosity) 
    << "set title 'Integral luminosity estimation ("+to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs)'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Integral luminosity, nb^{-1}'" 
    << "set xrange [0:45]"<< "set yrange [0:]";
}

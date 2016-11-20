// this file is distributed under 
// MIT license
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Genetic/searchmin.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
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
    vector<string> reaction={"He3eta","He3pi0pi0","He3pi0pi0pi0","He3pi0"};
    vector<hist<double>> norm;
    for(const string& r:reaction)norm.push_back(Hist(MC,r,histpath_forward_reconstr,"0-Reference"));

    vector<hist<double>> acceptance;
    for(const auto&h:norm)acceptance.push_back(hist<double>());
    hist<double> luminosity,bg_chi_sq,bg_ratio;
    RANDOM r_eng;
    for(size_t bin_num=0,bin_count=norm[0].size();bin_num<bin_count;bin_num++)
	if(norm[0][bin_num].X()>5.0){
	    auto Q=norm[0][bin_num].X();
	    string Qmsg="Q in ["+to_string(norm[0][bin_num].X().min())+":"+to_string(norm[0][bin_num].X().max())+"] MeV";
	    auto transform=[](hist<double>&h){h=h.XRange(0.45,0.58);};

	    hist<double> data=Hist(DATA,"",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
	    transform(data);
	    Plot<double> exp_plot;
	    exp_plot.Hist(data,"DATA")
	    << "set key on"<< "set title '"+Qmsg+"'"
	    << "set xlabel 'Missing mass, GeV'"
	    << "set ylabel 'counts'"
	    << "set yrange [0:]";
	    vector<hist<double>> theory;{
		Plot<double> th_plot;
		th_plot<< "set yrange [0:]"<< "set key on"<< "set xlabel 'Missing mass, MeV'"
		<< "set title '"+Qmsg+"'"
		<< "set ylabel 'acceptance density, GeV^{-1}'";
		for(size_t i=0;i<reaction.size();i++){
		    hist<double> react_sim=Hist(MC,reaction[i],histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
		    transform(react_sim);
		    auto N=norm[i][bin_num].Y();
		    auto MN=value<double>::std_error(react_sim.TotalSum().val());
		    acceptance[i] << point<value<double>>(Q,MN/N);
		    theory.push_back(react_sim/N);
		    th_plot.Hist(react_sim/(N*react_sim[0].X().uncertainty()*0.002),reaction[i]);
		}
	    }

	    SearchMin<DifferentialMutations<ParabolicErrorEstimationFromChisq>>
	    fit([&theory,&data](const ParamSet&P){
		double res=0;
		for(size_t i=0,n=data.size();i<n;i++){
		    value<double> exp_p=data[i].Y(),the_p=0;
		    for(size_t j=0,n=theory.size()-1;j<n;j++)
			the_p+=theory[j][i].Y()*P[j];
		    res+=exp_p.NumCompare(the_p);
		}
		return res;
	    });
	    const auto ex=parEq(3,0.0001);
	    fit.SetUncertaintyCalcDeltas(ex)
	    .SetFilter(make_shared<Above>()<<0.0<<0.0<<0.0);
	    const auto&data_count=data.TotalSum().val();
	    fit.Init(300,
		 make_shared<GenerateUniform>()
		     <<make_pair(0.0,20.0*data_count)
		     <<make_pair(0.0,20.0*data_count)
		     <<make_pair(0.0,20.0*data_count),
		 r_eng
	    );
	    while(
		!fit.AbsoluteOptimalityExitCondition(0.00001)
		&&
		!fit.ParametersDispersionExitCondition(ex)
	    ){
		fit.Iterate(r_eng);
		cout<<fit.iteration_count()<<" iterations; "
		<<fit.Optimality()<<"<chi^2<"
		<<fit.Optimality(fit.PopulationSize()-1)
		<<"        \r";
	    }
	    const auto&P=fit.ParametersWithUncertainties();
	    bg_ratio << point<value<double>>(Q,P[1]/P[2]);
	    bg_chi_sq << point<value<double>>(Q,fit.Optimality()/(data.size()-fit.ParamCount()));
	    exp_plot
	    .Line(hist<double>(theory[0]*P[0]+theory[1]*P[1]+theory[2]*P[2]).toLine(),"Total fit")
	    .Hist(theory[0]*P[0],"^3He eta")
	    .Hist(theory[1]*P[1],"^3He 3pi^0")
	    .Hist(theory[2]*P[2],"^3He 2pi^0");
		
	    luminosity << point<value<double>>(Q,
	       (P[0]/he3eta_sigma()(Q))
	       *double(trigger_he3_forward.scaling)
	    );
	}
    Plot<double>().Hist(bg_chi_sq) 
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'chi^2/d, n.d.'" 
    << "set yrange [0:]";

    {//Plot acceptance
	Plot<double> plot;
	plot << "set key on" 
	<< "set title 'Acceptance'"
	<< "set yrange [0:1.0]"
	<< "set xlabel 'Q, MeV'" 
	<< "set ylabel 'Acceptance, n.d.'";
	for(size_t i=0;i<reaction.size();i++)plot.Hist(acceptance[i],reaction[i]);
    }

    Plot<double>().Hist(bg_ratio) 
    << "set title 'N("+reaction[1]+")/N("+reaction[2]+")'"
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'n.d.'" 
    << "set yrange [0:]";

    auto runs=PresentRuns("");
    Plot<double>().Hist(luminosity) 
    << "set title 'Integral luminosity estimation ("+to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs)'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Integral luminosity, nb^{-1}'" 
    << "set yrange [0:]";

    Plot<double>()
    .Hist(hist<double>(he3eta_sigma().func(),BinsByStep(0.0,2.5,30.0)))
    .Hist(he3eta_sigma(),"Data from other experiments")
    << "set title 'Cross section of "+reaction[0]+" used in the calculations'"
    << "set key on" << "set xlabel 'Q, MeV'" 
    << "set ylabel 'sigma(^3He eta), nb'"
    << "set xrange [-10:50]"<< "set yrange [0:600]";
}

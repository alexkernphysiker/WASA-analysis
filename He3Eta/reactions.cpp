// this file is distributed under 
// GPL license
#include <iostream>
#include <iomanip>
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
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS),"background-reactions-forward");
    const auto runs=PresentRuns("F");
    const string runmsg=to_string(int(runs.first))+" of "+to_string(int(runs.second))+" runs";
    vector<string> histpath_forward_reconstr={"Histograms","He3Forward_Reconstruction"};
    vector<string> reaction={"He3eta","He3pi0pi0pi0","He3pi0pi0","He3pi0"};
    vector<hist<>> norm;
    for(const string& r:reaction)
	norm.push_back(Hist(MC,r,histpath_forward_reconstr,"0-Reference"));

    const hist<> luminosity=Plotter::Instance().GetPoints4("LUMINOSITYc");
    vector<hist<>> acceptance,true_events;
    for(size_t i=0;i<norm.size();i++){
	acceptance.push_back(hist<>());
	true_events.push_back(hist<>());
    }
    hist<> bg_chi_sq;
    vector<hist<>> fit_params;
    for(size_t i=0;i<reaction.size();i++)fit_params.push_back(hist<>());
    RANDOM r_eng;
    for(size_t bin_num=0,bin_count=norm[0].size();bin_num<bin_count;bin_num++){
	const auto&Q=norm[0][bin_num].X();
	const string Qmsg=static_cast<stringstream&>(stringstream()
	    <<"Q in ["<<setprecision(3)
	    <<Q.min()<<"; "<<Q.max()<<"] MeV"
	).str();
	auto transform=[](hist<>&h){h=h.XRange(0.40,0.60);};
	auto transform2=[](hist<>&h){h=h.XExclude(0.53,0.555);};
	hist<> data=Hist(DATA,"F",histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
	transform(data);
	Plot<>(Q.Contains(21)?"He3forward-above-mm":(Q.Contains(-39)?"He3forward-below-mm":"")).Hist(data)
	<< "set key on"<< "set title '"+Qmsg+", "+runmsg+"'"
	<< "set xlabel 'Missing mass, GeV'"
	<< "set ylabel 'counts'"
	<< "set yrange [0:]"<<"unset log y";
	Plot<> exp_plot(Q.Contains(21)?"He3forward-above-fit":(Q.Contains(-39)?"He3forward-below-fit":""));
	transform2(data);
	exp_plot.Hist(data,"DATA")
	<< "set key on"<< "set title '"+Qmsg+", "+runmsg+"'"
	<< "set xlabel 'Missing mass, GeV'"
	<< "set ylabel 'counts'"
	<< "set yrange [0:1500]"<<"unset log y";
	vector<hist<>> theory;{
	    Plot<> th_plot(Q.Contains(21)?"He3forward-above-mc":(Q.Contains(-39)?"He3forward-below-mc":""));
	    th_plot<< "set key on"<< "set xlabel 'Missing mass, MeV'"
	    << "set title '"+Qmsg+"'"
	    << "set ylabel 'acceptance density, GeV^{-1}'"
	    << "set yrange [0:20]";
	    for(size_t i=0;i<reaction.size();i++){
		hist<> react_sim=Hist(MC,reaction[i],histpath_forward_reconstr,string("MissingMass-Bin-")+to_string(bin_num));
		transform(react_sim);
		auto N=norm[i][bin_num].Y();
		th_plot.Line(react_sim.toLine()/(N.val()*react_sim[0].X().uncertainty()*2.),reaction[i]);
		auto MN=value<>::std_error(react_sim.TotalSum().val());
		transform2(react_sim);
		if(N.Above(0)){
		    acceptance[i] << point<value<>>(Q,MN/N);
		    theory.push_back(react_sim/N);
		}else{
		    acceptance[i] << point<value<>>(Q,0.0);
		    theory.push_back(react_sim*0.0);
		}
	    }
	}
	cout<<endl<<Qmsg<<endl<<endl;
	SearchMin<DifferentialMutations<Uncertainty>>
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
	fit.SetFilter([](const ParamSet&P)->bool{return (P[0]>=0)&&(P[1]>=0)&&(P[2]>=0)&&(P[3]>=0);});
	const auto&data_count=data.TotalSum().val();
	fit.Init(100,
	    make_shared<InitialDistributions>()
		<<make_shared<FixParam>(0)
		<<make_shared<DistribUniform>(0.0,0.1*data_count)
		<<make_shared<DistribUniform>(0.0,0.1*data_count)
		<<make_shared<FixParam>(0)
	    ,r_eng
	);
	while(!fit.AbsoluteOptimalityExitCondition(0.000000000001)){
	    fit.Iterate(r_eng);
	    cout<<fit.iteration_count()<<" iterations; "
	    <<fit.Optimality()<<"<chi^2<"
	    <<fit.Optimality(fit.PopulationSize()-1)
	    <<"          \r";
	}
	fit.SetUncertaintyCalcDeltas({0.1,0.1,0.1,0.1});
	const auto&P=fit.ParametersWithUncertainties();
	for(size_t i=0;i<reaction.size();i++)
	    fit_params[i]<< point<value<double>>(Q,P[i]);
	bg_chi_sq << point<value<double>>(Q,fit.Optimality()/(data.size()-2));
	exp_plot
	    .Hist(theory[1]*P[1]+theory[2]*P[2],"Background")
	    .Hist(theory[1]*P[1],reaction[1])
	    .Hist(theory[2]*P[2],reaction[2])
	;
	for(size_t i=0;i<reaction.size();i++)
	    true_events[i]<< point<value<double>>(Q,P[i]*double(trigger_he3_forward.scaling));
    }
    Plot<double>("He3forward-chisq").Hist(bg_chi_sq) 
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'chi^2/d, n.d.'" 
    << "set yrange [0:]"<<"unset log y";

    Plot<double> acc("He3forward-acceptance"),par,cs("He3forward-bg-cs");
    acc << "set key on"
    << "set yrange [0:1.0]"<<"unset log y"
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'Acceptance, n.d.'";
    par << "set key on" 
    << "set title 'Fitted coefficients'"
    << "set yrange [0:]"
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'coefficient, n.d.'";
    cs << "set key on" 
    << "set yrange [0:]"
    << "set xlabel 'Q, MeV'" 
    << "set ylabel 'cross section, nb'";
    for(size_t i=0;i<reaction.size();i++){
	acc.Hist(acceptance[i],reaction[i]);
	par.Hist(fit_params[i],reaction[i]);
	if((i>0)&&(i<3))cs.Hist(true_events[i]*runs.second/runs.first/luminosity,reaction[i],"CS-"+reaction[i]);
    }
}

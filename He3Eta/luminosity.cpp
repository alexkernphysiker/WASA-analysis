// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <math_h/sigma3.h>
#include <Genetic/fit.h>
#include <Genetic/uncertainties.h>
#include <Genetic/paramfunc.h>
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
#include <Parameters/parameters.h>
#include <Parameters/systematic.h>
#include "he3eta.h"
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "luminosity-forward");
    vector<string> histpath_forward_reconstr = {"Histograms", "He3Forward_Reconstruction"};
    const auto runs = PresentRuns("All");
    const string runmsg = (runs.first==runs.second)?"":"("+to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs)";
    const hist<> norm = Hist(MC, "He3eta-gg", histpath_forward_reconstr, "0-Reference");
    ext_hist<2> events_count,acceptance;
    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++)
        if (norm[bin_num].X() > 10.) {
            const auto &Q = norm[bin_num].X();
            const auto &N = norm[bin_num].Y();
            const string Qmsg =
                static_cast<stringstream &>(stringstream()
                    << "Q_{3Heη} є [" << setprecision(3)
                    << Q.min() << "; " << Q.max() << "] MeV"
                ).str();
            const string hist_name=string("MissingMass-Bin-") + to_string(bin_num);
            const list<size_t> params{pbeam_corr,he3_cut_h,he3_theta_cut};
            acceptance << make_point(Q,RawSystematicError(params,[&N,&histpath_forward_reconstr,&hist_name]
                (const string&suffix){
                const auto mc_unnorm = Hist(MC, "He3eta-gg", histpath_forward_reconstr,hist_name,suffix).XRange(0.530, 0.559);
                const auto mc = mc_unnorm / N;
                return extend_value<2,2>(mc.TotalSum());
            })());
            size_t counter=0;
            events_count << make_point(Q,RawSystematicError(params,[&counter,&Q,&N,&Qmsg,&histpath_forward_reconstr,&hist_name,&runmsg]
                (const string&suffix){
                counter++;
                const auto data_full = Hist(DATA, "All", histpath_forward_reconstr,hist_name,suffix).XRange(0.535, 0.559);
                const double bg_level=data_full.TransponateAndSort().right().X().max()*0.04;
                const auto data = data_full.XRange(0.535, data_full.YRange(bg_level,INFINITY).right().X().val()+0.001);
                const auto mc_unnorm = Hist(MC, "He3eta-gg", histpath_forward_reconstr,hist_name,suffix).XRange(0.535, 0.559);
                const auto chain = ChainWithStep(0.535, 0.001, 0.559);
                const auto mc = mc_unnorm / N;
                if(counter==1)Plot(Q.Contains(21) ? "He3eta-mc" : "",5)
                    .Hist(mc)
                        << "set key on" << "set title '" + Qmsg + ",3He+eta MC'"
                        << "set xlabel '3He missing mass, GeV'"
                        << "set ylabel 'Efficiency density, GeV^{-1}'"
                        << "set yrange [0:]" << "unset log y";
                cout << endl << Qmsg<< suffix << endl << endl;
                function<Uncertainties<2>(const double&,const double&)> func=
                                    [&data_full,&data,&counter,&Q,&Qmsg,&runmsg,&chain,&mc]
                    (const double&x,const double&y){
                    counter++;
                    const auto &data_count = data.TotalSum().val();
                    const auto data_bg = data.XExclude(x,y);
                    Fit<DifferentialMutations<>,ChiSquare,FunctionUncertaintiesEstimation> FIT(
                        removeXerorbars(data_bg),
                        [&data_count](const ParamSet & X, const ParamSet & P) {
                            return data_count * Polynom<4>(X[0], P);
                        }
                    );
                    FIT.SetUncertaintyCalcDeltas({0.1, 0.1,0.1,0.01,0.01})
                    .SetFilter([&FIT,&data](const ParamSet&P){
                        return
                        (FIT.func({data.left().X().min()},P)>=data.left().Y().min())&&
                        (FIT.func({data.right().X().val()},P)<=data.right().Y().max())&&
                        (FIT.func({data.right().X().min()},P)>0);
                    });
                    FIT.Init(300,
                        make_shared<InitialDistributions>()
                            << make_shared<DistribGauss>(0, 20)
                            << make_shared<DistribGauss>(0, 10)
                            << make_shared<DistribGauss>(0, 5)
                            << make_shared<DistribGauss>(0, 5)
                            << make_shared<DistribGauss>(0, 1)
                    );
                    cout << endl;
                    while (!FIT.AbsoluteOptimalityExitCondition(0.00001))FIT.Iterate();
                    cout << "Fitting: " << FIT.iteration_count() << " iterations; "
                        << FIT.Optimality() << "<chi^2<"
                        << FIT.Optimality(FIT.PopulationSize() - 1)
                        << endl;
                    const auto &P = FIT.ParametersWithUncertainties();
                    cout << endl;
                    const hist<> bg([&FIT](const value<>&X){
                        const auto r=FIT.FuncWithUncertainties({X.val()});
                        return (r>0)?r:value<>(0);
                    },data_full);
                    const auto clean = data_full - bg,clean2=clean.XRange(x,y);
                    if(counter==2){
                        const auto max=to_string(data_full.TransponateAndSort().right().X().max()*1.5/1000.);
                        Plot exp_plot(Q.Contains(21) ? "He3eta-fit" : (Q.Contains(14) ? "He3eta-fit-lo":""),5);
                        exp_plot.Hist(data_full/1000.,"Data").Hist(data_bg/1000.)
                            << "set key on" << "set title '" + Qmsg + runmsg + "'"
                            << "set xlabel '^3He missing mass, GeV/c^2'"
                            << "set ylabel 'Events count, 10^3'"
                            << "set yrange [-0.2:"+max+"]" << "unset log y";
                        const SortedPoints<> background([&FIT](double x) {return FIT({x});}, chain);
                        exp_plot.Line(background.YRange(0,INFINITY)/1000.).Hist(bg/1000.,"BG");
                        Plot subplot(Q.Contains(21) ? "He3eta-subtract" : (Q.Contains(14) ? "He3eta-subtract-lo":""),5);
                        subplot.Hist(clean/1000.).Line(Points<>{{clean.left().X().min(), 0.0},{clean.right().X().max(), 0.0}});
                        subplot.Hist(clean2/1000.,"DATA-BG")
                            .Line(toLine(mc*clean2.TotalSum()/mc.TotalSum())/1000.,"MC")
                            << "set key on" << "set title '" + Qmsg + runmsg + "'"
                            << "set xlabel '3He missing mass, GeV/c^2'"
			    << "set ylabel 'Events count, 10^3'"
                            << "set yrange [-0.2:"+max+"]" << "unset log y";
                    }
                    return extend_value<1,2>(clean.TotalSum()) / extend_value<2,2>(mc.TotalSum());
                };
                return SystematicError<he3eta_cut_left,he3eta_cut_right>(func)();
            })());
        }

    const auto cross_section=hist<>(he3eta_sigma().func(), BinsByStep(2.5, 2.5, 30.0));
    const hist<> exp_data=he3eta_sigma().XRange(11,35);
    Plot("He3eta-cross-section",5)
    .Hist(exp_data).Points(toLine(exp_data), "Experimental data","","with points pointtype 7 pointsize 3")
    .Hist(cross_section.XRange(10,35), "Interpolation", "CS-He3eta-assumed")
            << "set title 'Cross section'"
            << "set key on" << "set xlabel 'Q_{3Heη}, MeV'"
            << "set ylabel 'σ, nb'"
            << "set xrange [10:35]" << "set yrange [0:600]";

    Plot("He3eta-acceptance",5)
        .Hist(wrap_hist(acceptance.XRange(10,30)),"")
            << "set key on" << "set xlabel 'Q_{3Heη}, MeV'"
            << "set ylabel 'Efficiency, n.d.'"
            << "set xrange [10:30]" << "set yrange [0:1]";

    const auto luminosity=events_count*trigger_he3_forward.scaling
		/extend_hist<2,2>(cross_section.XRange(events_count.left().X().min(),events_count.right().X().max()));
    const auto Luminosity=add_one_uncertainty(events_count)*trigger_he3_forward.scaling
		/extend_hist<3,3>(cross_section.XRange(events_count.left().X().min(),events_count.right().X().max()));
    Plot("He3eta-luminosity",5)
        .Hist_2bars<1,2>(luminosity.XRange(12,30),"stat","syst","LUMINOSITYf")
            << "set title 'Integrated luminosity " + runmsg + "'"
            << "set key on" << "set xlabel 'Q_{3Heη}, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [10:30]" << "set yrange [0:]";
    cout<<"Partial estimation:"<<Luminosity.XRange(12.5,30).TotalSum([](Uncertainties<3>&s){
		s.use_maximum_estimation<2>();
		s.use_maximum_estimation<3>();
	})<<endl;
}

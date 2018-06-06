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
#include <Genetic/paramfunc.h>
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
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "luminosity-forward");
    vector<string> histpath_forward_reconstr = {"Histograms", "He3Forward_Reconstruction"};
    const auto runs = PresentRuns("All");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    const hist<> norm = Hist(MC, "He3eta-gg", histpath_forward_reconstr, "0-Reference");
    hist<> data_chi_sq, acceptance;
    ext_hist<2> events_count;
    vector<hist<>> parhists;
    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++)
        if (norm[bin_num].X() > 10.0) {
            const auto &Q = norm[bin_num].X();
            const auto &N = norm[bin_num].Y();
            const string Qmsg =
                static_cast<stringstream &>(stringstream()
                    << "Q in [" << setprecision(3)
                    << Q.min() << "; " << Q.max() << "] MeV"
                ).str();
            const string hist_name=string("MissingMass-Bin-") + to_string(bin_num);
            const auto data_full = Hist(DATA, "All", histpath_forward_reconstr,hist_name).XRange(0.532, 0.565);
            const double bg_level=data_full.TransponateAndSort().right().X().max()*0.04;
            const auto data = data_full.XRange(0.532, data_full.YRange(bg_level,INFINITY).right().X().val()+0.001);
            const auto mc_unnorm = Hist(MC, "He3eta-gg", histpath_forward_reconstr,hist_name).XRange(0.532, 0.565);
            const auto chain = ChainWithStep(0.532, 0.001, 0.565);
            const auto cut = make_pair(0.543,0.556);
            const auto mc = mc_unnorm / N;
            acceptance << make_point(Q, mc.TotalSum());
            Plot(Q.Contains(21) ? "He3eta-mc" : "")
            .Hist(mc)
                    << "set key on" << "set title '" + Qmsg + ",3He+eta MC'"
                    << "set xlabel 'Missing mass, GeV'"
                    << "set ylabel 'acceptance density, GeV^{-1}'"
                    << "set yrange [0:]" << "unset log y";
            const auto &data_count = data.TotalSum().val();
            const auto data_bg = data.XExclude(cut.first, cut.second);
            cout << endl << Qmsg << endl << endl;
            Fit2<DifferentialMutations<>> FIT(
                data_bg.removeXerorbars(),
                [&data_count](const ParamSet & X, const ParamSet & P) {
                    return data_count * Polynom<4>(X[0], P);
                }
            );
            FIT.SetUncertaintyCalcDeltas({0.1, 0.1, 0.1, 0.1, 0.01})
            .SetFilter([&FIT,&data,&Q](const ParamSet&P){
                return
                (FIT.func({data.right().X().val()},P)<=data.right().Y().max())&&
                (FIT.func({data.left().X().min()},P)>=data.left().Y().min())&&
                ((Q<10.)||(FIT.func({0.555},P)>0));
            });
            FIT.Init(200,
                     make_shared<InitialDistributions>()
                        << make_shared<DistribGauss>(0, 10)
                        << make_shared<DistribGauss>(0, 10)
                        << make_shared<DistribGauss>(0, 5)
                        << make_shared<DistribGauss>(0, 5)
                        << make_shared<DistribGauss>(0, 1)
            );
            cout << endl;
            while (!FIT.AbsoluteOptimalityExitCondition(0.000001))FIT.Iterate();
            cout << "Fitting: " << FIT.iteration_count() << " iterations; "
                 << FIT.Optimality() << "<chi^2<"
                 << FIT.Optimality(FIT.PopulationSize() - 1)
                 << endl;
            const auto &P = FIT.ParametersWithUncertainties();
            if (parhists.size() == 0) {
                for (size_t i = 0; i < P.size(); i++)
                    parhists.push_back(hist<>());
            }
            for (size_t i = 0; i < P.size(); i++)
                parhists[i] << make_point(Q, P[i]);
            data_chi_sq << make_point(Q, FIT.Optimality() / (data.size() - FIT.ParamCount()));
            cout << endl;
            Plot exp_plot(Q.Contains(21) ? "He3eta-fit" : (Q.Contains(14) ? "He3eta-fit-lo":""));
            exp_plot.Hist(data_full).Hist(data_bg)
                    << "set key on" << "set title '" + Qmsg + ", " + runmsg + "'"
                    << "set xlabel 'Missing mass, GeV'"
                    << "set ylabel 'counts'"
                    << "set yrange [-200:]" << "unset log y";
            const SortedPoints<> background([&FIT](double x) {return FIT({x});}, chain);
            const hist<> bg([&FIT](const value<>&X){
                const auto r=FIT.FuncWithUncertainties({X.val()});
                return (r>0)?r:value<>(0);
            },data_full);
            exp_plot.Line(background.YRange(0,INFINITY)).Hist(bg);
            auto clean = data_full - bg;
            Plot subplot(Q.Contains(21) ? "He3eta-subtract" : (Q.Contains(14) ? "He3eta-subtract-lo":""));
            subplot.Hist(clean,"DATA").Line(Points<>{{clean.left().X().min(), 0.0},{clean.right().X().max(), 0.0}});
            subplot.Hist(clean = clean.XRange(cut.first, cut.second))
                .Line((mc*clean.TotalSum()/mc.TotalSum()).toLine(),"MC")
                    << "set key on" << "set title '" + Qmsg + ", " + runmsg + "'"
                    << "set xlabel 'Missing mass, GeV'"
                    << "set ylabel 'counts'"
                    << "set yrange [-200:]" << "unset log y";

            events_count << make_point(Q,extend_value<1,2>(clean.TotalSum()) / extend_value<2,2>(mc.TotalSum()));
        }
    for (size_t i = 0; i < parhists.size(); i++)
        Plot().Hist(parhists[i])
                << "set xlabel 'Q, MeV'"
                << "set ylabel 'parameter" + to_string(i) + "'";
    Plot("He3eta-chisq").Hist(data_chi_sq)
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'chi^2/d, n.d.'"
            << "set yrange [0:]" << "unset log y";

    const auto cross_section=hist<>(he3eta_sigma().func(), BinsByStep(2.5, 2.5, 30.0));
    Plot("He3eta-cross-section")
    .Hist(he3eta_sigma().XRange(12.5,35), "Experimental data")
    .Hist(cross_section.XRange(12.5,35), "Interpolation", "CS-He3eta-assumed")
            << "set title 'Cross section of He3eta used in the calculations'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'sigma(^3He eta), nb'"
            << "set xrange [10:35]" << "set yrange [0:600]";

    Plot("He3eta-acceptance")
        .Hist(acceptance.XRange(12.5,30),"")
            << "set title '3He+eta acceptance'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'acceptance, n.d.'"
            << "set xrange [10:30]" << "set yrange [0:1]";

    const auto luminosity=events_count*trigger_he3_forward.scaling/extend_hist<2,2>(cross_section.XRange(events_count.left().X().min(),events_count.right().X().max()));
    Plot("He3eta-luminosity")
        .Hist_2bars<1,2>(luminosity.XRange(12.5,30),"stat","syst","LUMINOSITYf")
            << "set title 'Integrated luminosity (" + runmsg + ")'"
            << "set key on" << "set xlabel 'Q, MeV'"
            << "set ylabel 'Integrated luminosity, nb^{-1}'"
            << "set xrange [10:30]" << "set yrange [0:]";
}

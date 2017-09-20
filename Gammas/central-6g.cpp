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
#include <Genetic/initialconditions.h>
#include <Genetic/filter.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
using namespace std;
using namespace ROOT_data;
using namespace Genetic;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    Plotter<>::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-6gamma-with-3he");
    vector<string> histpath_reconstr = {"Histograms", "He3nCentralGammas"};
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas6"};
    vector<string> reaction = {"bound1-6g", "He3eta", "He3pi0pi0pi0"};
    hist<> norm = Hist(MC, reaction[0], histpath_reconstr, "0-Reference");
    const auto runs = PresentRuns("C");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    Plot<> theory("He3gggggg-IMDiff0-mc"), theory1("He3gggggg-IMDiff1-mc"),
         experiment("He3gggggg-IMDiff-data");
    for (const auto &r : reaction) {
        theory.Line(Hist(MC, r, histpath_central_reconstr, "GIMDiff2").toLine(), r);
        theory1.Line(Hist(MC, r, histpath_central_reconstr, "GIMDiff3").toLine(), r);
    }
    theory << "set key on" << "set yrange [0:]";
    theory1 << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GIMDiff2"), "3He 6gamma")
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GIMDiff3"), "3He 3pi^0")
            << "set key on" << "set title '" + runmsg + "'" << "set yrange [0:]";
    Plot<> mm_theory("He3gggggg-MM0-mc"), mm_theory1("He3gggggg-MM1-mc")
    , mm_experiment("He3gggggg-MM-data");
    for (const auto &r : reaction) {
        mm_theory.Line(Hist(MC, r, histpath_central_reconstr, "GMM2").toLine(), r);
        mm_theory1.Line(Hist(MC, r, histpath_central_reconstr, "GMM3").toLine(), r);
    }
    mm_theory << "set key on" << "set yrange [0:]";
    mm_theory1 << "set key on" << "set yrange [0:]";
    mm_experiment
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GMM2"), "3He 6gamma")
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GMM3"), "3He 3pi^0")
            << "set key on" << "set title '" + runmsg + "'" << "set yrange [0:]";

    hist<> ev_am;
    vector<hist<>> acceptance;
    for (size_t i = 0; i < reaction.size(); i++) {
        acceptance.push_back(hist<>());
    }
    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
                            << "Q in [" << setprecision(3)
                            << Q.min() << "; " << Q.max() << "] MeV").str();
        Plot<> mc_plot(
            Q.Contains(21) ? "He3gggggg-above-mc" : (
                Q.Contains(-39) ? "He3gggggg-below-mc" : (
                    Q.Contains(-3) ? "He3gggggg-thr-mc" : ""
                )
            )
        );
        mc_plot << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]"
                << "set xlabel 'pi0+pi0+pi0 invariant mass, GeV'";
        Plot<> mc_plot1(
            Q.Contains(21) ? "He3gggggg-above-mc" : (
                Q.Contains(-39) ? "He3gggggg-below-mc" : (
                    Q.Contains(-3) ? "He3gggggg-thr-mc" : ""
                )
            )
        );
        mc_plot1 << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]"
                 << "set xlabel 'pi0+pi0+pi0 invariant mass, GeV'";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            hist<double> Norm = Hist(MC, r, histpath_reconstr, "0-Reference");
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) {
                const hist<> h = Hist(MC, r, histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num));
                const hist<> h1 = Hist(MC, r, histpath_central_reconstr, string("GIM3-Bin-") + to_string(bin_num));
                const auto C = h1.TotalSum();
                mc_plot.Hist(h / N, r);
                mc_plot1.Hist(h1 / N, r);
                acceptance[i] << point<value<double>>(Q, C / N);
            } else {
                acceptance[i] << point<value<>>(Q, 0.0);
            }
        }
        const hist<> data = Hist(DATA, "C", histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num));
        const hist<> data1 = Hist(DATA, "C", histpath_central_reconstr, string("GIM3-Bin-") + to_string(bin_num));
        Plot<>(
            Q.Contains(21) ? "He3gggggg-above-data" : (
                Q.Contains(-39) ? "He3gggggg-below-data" : (
                    Q.Contains(-3) ? "He3gggggg-thr-data" : ""
                )
            )
        )
        .Hist(data, "3He 6gamma") .Hist(data1, "3He 3pi^0") //.Hist(data2,"3pi^0 cut")
                << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                << "set xlabel 'pi0+pi0+pi0 invariant mass, GeV'";
        ev_am << point<value<>>(Q, value<>::std_error(data1.TotalSum().val()));
    }
    Plot<> accplot("He3gggggg-acceptance");
    accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:]" << "set key on";
    for (size_t i = 0; i < reaction.size(); i++) {
        const auto acc = acceptance[i].YRange(0.0001, INFINITY);
        if (acc.size() > 0) {
            accplot.Hist(acc, reaction[i]);
        }
    }
    const hist<> luminosity = Plotter<>::Instance().GetPoints<value<>>("LUMINOSITYc");
    const hist<> he3etacs = Plotter<>::Instance().GetPoints<value<>>("CS-He3eta-assumed");
    const hist<> known_events =
        luminosity * (runs.first / runs.second) / double(trigger_he3_forward.scaling)
        * (he3etacs * acceptance[1]);
    Plot<>("He3gggggg-events")
    .Hist(ev_am - known_events, "data-(3He+eta)")
    .Hist(ev_am, "data")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'";
}

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
    Plotter<>::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-6gamma-only");
    vector<string> histpath_reconstr = {"Histograms", "OnlyCentralGammas"};
    vector<string> histpath_central_reconstr = {"Histograms", "OnlyCentralGammas6"};
    vector<string> reaction = {"He3eta", "He3pi0pi0pi0"};
    const hist<> norm = Hist(MC, reaction[0], histpath_reconstr, "0-Reference");
    const auto runs = PresentRuns("CC");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    Plot<> theory("gggggg-IMDiff0-mc"), theory1("gggggg-IMDiff1-mc"), theory2("gggggg-IMDiff2-mc"), experiment("gggggg-IMDiff-data");
    for (const auto &r : reaction) {
        theory.Line(Hist(MC, r, histpath_central_reconstr, "GIMPDiff0").toLine(), r);
        theory1.Line(Hist(MC, r, histpath_central_reconstr, "GIMPDiff1").toLine(), r);
        theory2.Line(Hist(MC, r, histpath_central_reconstr, "GIMPDiff2").toLine(), r);
    }
    theory << "set key on" << "set yrange [0:]";
    theory1 << "set key on" << "set yrange [0:]";
    theory2 << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "CC", histpath_central_reconstr, "GIMPDiff0"))
    .Hist(Hist(DATA, "CC", histpath_central_reconstr, "GIMPDiff1"))
    .Hist(Hist(DATA, "CC", histpath_central_reconstr, "GIMPDiff2"))
            << "set key on" << "set title '" + runmsg + "'" << "set yrange [0:]";
    Plot<> mm_theory("gggggg-MM0-mc"), mm_theory1("gggggg-MM1-mc"), mm_theory2("gggggg-MM1-mc"), mm_experiment("gggggg-MM-data");
    for (const auto &r : reaction) {
        mm_theory.Line(Hist(MC, r, histpath_central_reconstr, "GMM0").toLine(), r);
        mm_theory1.Line(Hist(MC, r, histpath_central_reconstr, "GMM1").toLine(), r);
        mm_theory2.Line(Hist(MC, r, histpath_central_reconstr, "GMM2").toLine(), r);
    }
    mm_theory << "set key on" << "set yrange [0:]";
    mm_theory1 << "set key on" << "set yrange [0:]";
    mm_theory2 << "set key on" << "set yrange [0:]";
    mm_experiment
    .Hist(Hist(DATA, "CC", histpath_central_reconstr, "GMM0"), "6gamma")
    .Hist(Hist(DATA, "CC", histpath_central_reconstr, "GMM1"), "3pi^0")
    .Hist(Hist(DATA, "CC", histpath_central_reconstr, "GMM2"), "MM cut")
            << "set key on" << "set title '" + runmsg + "'" << "set yrange [0:]";
    hist<> ev_am;
    vector<hist<>> acceptance;
    for (size_t i = 0; i < reaction.size(); i++) {
        acceptance.push_back(hist<>());
    }
    for (size_t bin_num = 0; bin_num < norm.size(); bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
                            << "Q in [" << setprecision(3)
                            << Q.min() << "; " << Q.max() << "] MeV").str();
        Plot<> mc_plot(
            Q.Contains(21) ? "gggggg-above-mc0" : (
                Q.Contains(-39) ? "gggggg-below-mc0" : (
                    Q.Contains(-3) ? "gggggg-thr-mc0" : ""
                )
            )
        );
        Plot<> mc_plot1(
            Q.Contains(21) ? "gggggg-above-mc1" : (
                Q.Contains(-39) ? "gggggg-below-mc1" : (
                    Q.Contains(-3) ? "gggggg-thr-mc1" : ""
                )
            )
        );
        Plot<> mc_plot2(
            Q.Contains(21) ? "gggggg-above-mc2" : (
                Q.Contains(-39) ? "gggggg-below-mc2" : (
                    Q.Contains(-3) ? "gggggg-thr-mc2" : ""
                )
            )
        );
        mc_plot << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]";
        mc_plot1 << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]";
        mc_plot2 << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            hist<> Norm = Hist(MC, r, histpath_reconstr, "0-Reference");
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) {
                const hist<> h = Hist(MC, r, histpath_central_reconstr, string("GIM0-Bin-") + to_string(bin_num));
                const hist<> h1 = Hist(MC, r, histpath_central_reconstr, string("GIM1-Bin-") + to_string(bin_num));
                const hist<> h2 = Hist(MC, r, histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num));
                const auto C = h2.TotalSum();
                mc_plot.Hist(h / N, r);
                mc_plot1.Hist(h1 / N, r);
                mc_plot2.Hist(h2 / N, r);
                acceptance[i] << point<value<>>(Q, C / N);
            } else {
                acceptance[i] << point<value<>>(Q, 0.0);
            }
        }
        const hist<> data = Hist(DATA, "CC", histpath_central_reconstr, string("GIM0-Bin-") + to_string(bin_num));
        const hist<> data1 = Hist(DATA, "CC", histpath_central_reconstr, string("GIM1-Bin-") + to_string(bin_num));
        const hist<> data2 = Hist(DATA, "CC", histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num));
        Plot<>(
            Q.Contains(21) ? "gggggg-above-data" : (
                Q.Contains(-39) ? "gggggg-below-data" : (
                    Q.Contains(-3) ? "gggggg-thr-data" : ""
                )
            )
        )
        .Hist(data, "All").Hist(data1, "3pi^0").Hist(data2, "MM cut")
                << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]" << "set key on";
        ev_am << point<value<>>(Q, value<>::std_error(data2.TotalSum().val()));
    }
    Plot<> accplot("gggggg-acceptance");
    accplot
            << "set title 'Acceptance'"
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
    const hist<> he3eta_events = luminosity * (runs.first / runs.second) 
            / double(trigger_gammas_central.scaling) 
            * (he3etacs * acceptance[0]);
    const hist<> he3pi0pi0pi0_events = luminosity * (runs.first / runs.second) 
            / double(trigger_gammas_central.scaling) 
            * (acceptance[1]*value<>(115,25));
    Plot<>("gggggg-events")
    .Hist(ev_am , "data")
    .Line(he3eta_events.toLine(), "3He+eta")
    .Line(he3pi0pi0pi0_events.toLine(), "3He+3pi0")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'"
            << "set title 'Events 6gamma'" << "set yrange [0:]";
    const hist<> cross_section = (ev_am - he3eta_events) / acceptance[1]
                                 * double(trigger_gammas_central.scaling)
                                 * (runs.second / runs.first) / luminosity;
}

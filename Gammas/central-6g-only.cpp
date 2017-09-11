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
    vector<string> reaction = {"bound1-6g", "He3eta", "He3pi0pi0pi0"};
    const hist<> norm = Hist(MC, reaction[0], histpath_reconstr, "0-Reference");
    const auto runs = PresentRuns("CC");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    Plot<> theory("gggggg-IMDiff0-mc"),theory1("gggggg-IMDiff1-mc"), experiment("gggggg-IMDiff-data");
    for (const auto &r : reaction) {
        theory.Line(Hist(MC, r, histpath_central_reconstr, "GMMPDiff0").toLine(), r);
        theory1.Line(Hist(MC, r, histpath_central_reconstr, "GMMPDiff1").toLine(), r);
    }
    theory << "set key on" << "set yrange [0:]";
    theory1 << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "CC", histpath_central_reconstr, "GMMPDiff0"))
    .Hist(Hist(DATA, "CC", histpath_central_reconstr, "GMMPDiff1"))
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
            Q.Contains(21) ? "gggggg-above-mc" : (
                Q.Contains(-39) ? "gggggg-below-mc" : (
                    Q.Contains(-3) ? "gggggg-thr-mc" : ""
                )
            )
        );
        mc_plot << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            hist<> Norm = Hist(MC, r, histpath_reconstr, "0-Reference");
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) {
                const hist<> h = Hist(MC, r, histpath_central_reconstr, string("GIM1-Bin-") + to_string(bin_num));
                const auto C = h.TotalSum();
                mc_plot.Hist(h / N, r);
                acceptance[i] << point<value<>>(Q, C / N);
            } else {
                acceptance[i] << point<value<>>(Q, 0.0);
            }
        }
        const hist<>
        data = Hist(DATA, "CC", histpath_central_reconstr, string("GIM1-Bin-") + to_string(bin_num));
        Plot<>(
            Q.Contains(21) ? "gggggg-above-data" : (
                Q.Contains(-39) ? "gggggg-below-data" : (
                    Q.Contains(-3) ? "gggggg-thr-data" : ""
                )
            )
        )
        .Hist(data) << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]";
        ev_am << point<value<>>(Q, value<>::std_error(data.TotalSum().val()));
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
    const auto luminosity = Plotter<>::Instance().GetPoints4("LUMINOSITYc");
    const auto he3etacs = Plotter<>::Instance().GetPoints4("CS-He3eta-assumed");
    const hist<> he3eta_events = luminosity * (runs.first / runs.second) / double(trigger_gammas_central.scaling) * (he3etacs * acceptance[1]);
    Plot<>("gggggg-events")
    .Hist(he3eta_events, "3He+eta")
    .Hist(ev_am , "data")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'"
            << "set title 'Events 6gamma'" << "set yrange [0:]";
    const hist<> cross_section = (ev_am - he3eta_events)/acceptance[2]
        * double(trigger_gammas_central.scaling) 
        * (runs.second / runs.first) / luminosity;
    Plot<>("gggggg-cross-section")
    .Hist(cross_section , "data")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'cross section, nb'"
            << "set title 'He3+3pi0'" << "set yrange [0:]";
}

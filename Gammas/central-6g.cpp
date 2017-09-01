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
    Plotter<>::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-6gamma");
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas6"};
    vector<string> reaction = {"bound1-6g", "He3eta", "He3pi0", "He3pi0pi0", "He3pi0pi0pi0"};
    hist<> norm = Hist(MC, reaction[0], {"Histograms", "He3nCentralGammas"}, "0-Reference");
    const auto runs = PresentRuns("C");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    Plot<> theory("He3gggggg-IMDiff-mc"), experiment("He3gggggg-IMDiff-data");
    for (const auto &r : reaction) {
        theory.Line(Hist(MC, r, histpath_central_reconstr, "GIMDiff1").toLine(), r);
    }
    theory << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GIMDiff1"), "DATA")
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
                            << Q.min() << "; " << Q.max() << "] MeV"
                                                       ).str();
        Plot<> mc_plot(
            Q.Contains(21) ? "He3gggggg-above-mc": (
                Q.Contains(-39) ? "He3gggggg-below-mc": (
                    Q.Contains(-3) ? "He3gggggg-thr-mc": ""
                )
            )
        );
        mc_plot << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            hist<double> Norm = Hist(MC, r, {"Histograms", "He3nCentralGammas"}, "0-Reference");
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) {
                const hist<> h = Hist(MC, r, histpath_central_reconstr, string("GIM1-Bin-") + to_string(bin_num));
                const auto C = h.TotalSum();
                mc_plot.Hist(h / N, r);
                acceptance[i] << point<value<double>>(Q, C / N);
            }
        }
        const hist<>
        data = Hist(DATA, "C", histpath_central_reconstr, string("GIM1-Bin-") + to_string(bin_num));
        Plot<>(
            Q.Contains(21) ? "He3gggggg-above-data" : (
                Q.Contains(-39) ? "He3gggggg-below-data" : (
                    Q.Contains(-3) ? "He3gggggg-thr-data" : ""
                )
            )
        )
        .Hist(data) << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]";
        ev_am << point<value<double>>(Q, data.TotalSum());
    }
    Plot<> accplot("He3gggggg-acceptance");
    accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:]" << "set key on";
    for (size_t i = 0; i < reaction.size(); i++) {
        if (acceptance[i].size() > 0) {
            accplot.Hist(acceptance[i], reaction[i]);
        }
    }
    Plot<>().Hist(ev_am) << "set title 'Data events'" << "set yrange [0:]";
    Plot<>("He3gggggg-cs").Hist((ev_am / acceptance[0])
                                *trigger_he3_forward.scaling)
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'cross section, nb'"
            << "set title 'Events norm'" << "set yrange [0:]";
}

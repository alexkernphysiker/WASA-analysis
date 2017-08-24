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
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-2gamma");
    vector<string> histpath_central_reconstr = {"Histograms", "CentralGammas2"};
    vector<string> reaction = {"bound1-2g", "He3eta", "He3pi0", "He3pi0pi0", "He3pi0pi0pi0"};
    const auto runs = PresentRuns("C");
    hist<> norm = Hist(MC, reaction[0], {"Histograms", "CentralGammas"}, "0-Reference");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    Plot<> theory("He3gg-mc"), experiment("He3gg-data");
    for (const auto &r : reaction) {
        theory.Line(Hist(MC, r, histpath_central_reconstr, "GIMDiff-AllBins").toLine(), r);
    }
    theory << "set key on" << "set yrange [0:]";
    experiment.Hist(Hist(DATA, "C", histpath_central_reconstr, "GIMDiff-AllBins"), "DATA")
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
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            hist<> Norm = Hist(MC, r, {"Histograms", "CentralGammas"}, "0-Reference");
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) {
                const hist<> h = Hist(MC, r, histpath_central_reconstr, string("GIMDiff-Bin-") + to_string(bin_num)),
                             H = h.XRange(0, 0.05);
                const auto C = H.TotalSum();
                if (C.Above(0)) {
                    Plot<>(
                        Q.Contains(21) ? "He3gg-above-mc" + r : (
                            Q.Contains(-39) ? "He3gg-below-mc" + r : (
                                Q.Contains(-3) ? "He3gg-thr-mc" + r : ""
                            )
                        )
                    )
                    .Hist(h).Hist(H) << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]";
                    acceptance[i] << point<value<double>>(Q, C / N);
                }
            }
        }
        const hist<> data = Hist(DATA, "C", histpath_central_reconstr, string("GIMDiff-Bin-") + to_string(bin_num)),
                     DATA = data.XRange(0, 0.05);
        Plot<>(
            Q.Contains(21) ? "He3gg-above-data" : (
                Q.Contains(-39) ? "He3gg-below-data" : (
                    Q.Contains(-3) ? "He3gg-thr-data" : ""
                )
            )
        )
        .Hist(data).Hist(DATA) << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]";
        ev_am << point<value<double>>(Q, DATA.TotalSum());
    }
    Plot<> accplot("He3gg-acceptance");
    accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:]" << "set key on";
    for (size_t i = 0; i < reaction.size(); i++) {
        if (acceptance[i].size() > 0) {
            accplot.Hist(acceptance[i], reaction[i]);
        }
    }
    const auto luminosity = Plotter::Instance().GetPoints4("LUMINOSITYc");
    Plot<>().Hist(ev_am) << "set title 'Data events'" << "set yrange [0:]";
    Plot<>("He3gg-cs").Hist((ev_am / acceptance[0] / luminosity)
                            *trigger_he3_forward.scaling, "", "CS-He3gg")
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'cross section, nb'"
            << "set title 'Events norm'" << "set yrange [0:]";
}



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
    Plotter<>::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-2gamma");
    vector<string> histpath_reconstr = {"Histograms", "He3nCentralGammas"};
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound1-2g", "He3eta", "He3pi0pi0", "bound1-6g", "He3pi0pi0pi0", "He3pi0"};
    const auto runs = PresentRuns("C");
    const hist<> norm = Hist(MC, reaction[0], histpath_reconstr, "0-Reference");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    Plot<> theory("He3gg-IMDiff-mc"), experiment("He3gg-IMDiff-data");
    for (const auto &r : reaction) {
        theory.Line(Hist(MC, r, histpath_central_reconstr, "GIMDiff3").toLine(), r);
    }
    theory << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GIMDiff2"), "3He 2gamma")
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GIMDiff3"), "2gamma IM cut")
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GIMDiff4"), "2gamma MM cut")
            << "set title '" + runmsg + "'" << "set yrange [0:]";
    Plot<> mm_theory("He3gg-MM0-mc"), mm_theory1("He3gg-MM1-mc"),
         mm_theory2("He3gg-MM2-mc"), mm_experiment("He3gg-MM-data");
    for (const auto &r : reaction) {
        if (r == reaction[0]) {
            Plot<>("He3gg-MM-bound-mc")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMM2"), "3He 2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMM3"), "2gamma IM cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMM4"), "2gamma MM cut")
                    << "set key on" << "set yrange [0:]";
        }
        mm_theory.Line(Hist(MC, r, histpath_central_reconstr, "GMM2").toLine(), r);
        mm_theory1.Line(Hist(MC, r, histpath_central_reconstr, "GMM3").toLine(), r);
        mm_theory2.Line(Hist(MC, r, histpath_central_reconstr, "GMM4").toLine(), r);
    }
    mm_theory << "set key on" << "set yrange [0:]";
    mm_theory1 << "set key on" << "set yrange [0:]";
    mm_theory2 << "set key on" << "set yrange [0:]";
    mm_experiment
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GMM2"), "3He 2gamma")
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GMM3"), "2gamma IM cut")
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GMM4"), "2 gamma MM cut")
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
        Plot<>(
            Q.Contains(21) ? "He3gg-above-he3mm-bound-mc" : (
                Q.Contains(-39) ? "He3gg-below-he3mm-bound-mc" : (
                    Q.Contains(-3) ? "He3gg-thr-he3mm-bound-mc" : ""
                )
            )
        )
        .Hist(Hist(MC, reaction[0], histpath_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
        .Hist(Hist(MC, reaction[0], histpath_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
        .Hist(Hist(MC, reaction[0], histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+2gamma")
        .Hist(Hist(MC, reaction[0], histpath_central_reconstr, string("He3MM3-Bin-") + to_string(bin_num)), "2gamma IM cut")
        .Hist(Hist(MC, reaction[0], histpath_central_reconstr, string("He3MM4-Bin-") + to_string(bin_num)), "2gamma MM cut")
                << "set key on" << "set title '" + Qmsg + ";MC " + reaction[0] + "'" << "set yrange [0:]"
                << "set xlabel '3He missing mass, GeV'";
        Plot<>(
            Q.Contains(21) ? "He3gg-above-he3mm-data" : (
                Q.Contains(-39) ? "He3gg-below-he3mm-data" : (
                    Q.Contains(-3) ? "He3gg-thr-he3mm-data" : ""
                )
            )
        )
        .Hist(Hist(DATA, "C", histpath_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
        .Hist(Hist(DATA, "C", histpath_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
        .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+2gamma")
        .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM3-Bin-") + to_string(bin_num)), "2gamma IM cut")
        .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM4-Bin-") + to_string(bin_num)), "2gamma MM cut")
                << "set key on" << "set title '" + Qmsg + ";MC " + reaction[0] + "'" << "set yrange [0:]"
                << "set xlabel '3He missing mass, GeV'";
        Plot<> he3_plot(
            Q.Contains(21) ? "He3gg-above-he3mm-mc" : (
                Q.Contains(-39) ? "He3gg-below-he3mm-mc" : (
                    Q.Contains(-3) ? "He3gg-thr-he3mm-mc" : ""
                )
            )
        );
        he3_plot << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]"
                 << "set xlabel '3He missing mass, GeV'";
        Plot<> mc_plot(
            Q.Contains(21) ? "He3gg-above-im-mc" : (
                Q.Contains(-39) ? "He3gg-below-im-mc" : (
                    Q.Contains(-3) ? "He3gg-thr-im-mc" : ""
                )
            )
        );
        mc_plot << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]"
                << "set xlabel 'gamma-gamma invariant mass, GeV'";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            hist<> Norm = Hist(MC, r, histpath_reconstr, "0-Reference");
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) {
                const hist<> h = Hist(MC, r, histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num));
                const hist<> he3 = Hist(MC, r, histpath_central_reconstr, string("He3MM4-Bin-") + to_string(bin_num));
                const auto C = value<>::std_error(h.TotalSum().val());
                mc_plot.Hist(h / N, r);
                he3_plot.Hist(he3 / N, r);
                acceptance[i] << point<value<>>(Q, C / N);
                Plot<>(
                    Q.Contains(21) ? "He3gg-above-ggim-mc" + r : (
                        Q.Contains(-39) ? "He3gg-below-bound-ggim-mc" + r : (
                            Q.Contains(-3) ? "He3gg-thr-ggim-mc" + r : ""
                        )
                    )
                )
                .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num)), "3He with 2 gammas")
                .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM3-Bin-") + to_string(bin_num)), "2 gamma IM cut")
                .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num)), "2 gamma MM cut")
                        << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                        << "set xlabel 'gamma-gamma invariant mass, GeV'";
            } else {
                acceptance[i] << point<value<>>(Q, 0.0);
            }
        }
        const hist<> data1 = Hist(DATA, "C", histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num));
        const hist<> data2 = Hist(DATA, "C", histpath_central_reconstr, string("GIM3-Bin-") + to_string(bin_num));
        const hist<> data3 = Hist(DATA, "C", histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num));
        Plot<>(
            Q.Contains(21) ? "He3gg-above-data" : (
                Q.Contains(-39) ? "He3gg-below-data" : (
                    Q.Contains(-3) ? "He3gg-thr-data" : ""
                )
            )
        )
        .Hist(data1, "All").Hist(data2, "3pi^0").Hist(data3, "MM cut")
                << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                << "set xlabel 'gamma-gamma invariant mass, GeV'";
        ev_am << point<value<>>(Q, value<>::std_error(data3.TotalSum().val()));
    }
    Plot<> accplot("He3gg-acceptance");
    accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:]" << "set key on";
    Plot<> accplot2("He3gg-acceptance2");
    accplot2 << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:0.001]" << "set key on";
    for (size_t i = 0; i < reaction.size(); i++) {
        const auto acc = acceptance[i].YRange(0.000, INFINITY);
        if (acc.size() > 0) {
            accplot.Hist(acc, reaction[i]);
            accplot2.Hist(acc, reaction[i]);
        }
    }
    const hist<> luminosity = Plotter<>::Instance().GetPoints<value<>>("LUMINOSITYc");
    const hist<> he3etacs = Plotter<>::Instance().GetPoints<value<>>("CS-He3eta-assumed");
    const hist<> known_events =
        luminosity * (runs.first / runs.second) / double(trigger_he3_forward.scaling)
        * (he3etacs * acceptance[1]);
    Plot<>("He3gg-events")
    .Hist(ev_am, "data")
    .Line(known_events.toLine(), "3He+eta estimated")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'";
}

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
    vector<string> reaction = {"bound1-2g","He3eta","He3pi0pi0", "bound1-6g","He3pi0pi0pi0", "He3pi0"};
    const auto runs = PresentRuns("C");
    const hist<> norm = Hist(MC, reaction[0], histpath_reconstr, "0-Reference");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    Plot<> theory("He3gg-IMDiff-mc"), experiment("He3gg-IMDiff-data");
    for (const auto &r : reaction) {
        theory.Line(Hist(MC, r, histpath_central_reconstr, "GIMDiff1").toLine(), r);
    }
    theory << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GIMDiff0"), "DATA before cut")
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GIMDiff1"), "DATA after cut")
            << "set title '" + runmsg + "'" << "set yrange [0:]";
    Plot<> mm_theory("He3gg-MM0-mc"),mm_theory1("He3gg-MM1-mc"), mm_experiment("He3gg-MM-data");
    for (const auto &r : reaction) {
        mm_theory.Line(Hist(MC, r, histpath_central_reconstr, "GMM0").toLine(), r);
        mm_theory1.Line(Hist(MC, r, histpath_central_reconstr, "GMM1").toLine(), r);
    }
    mm_theory << "set key on" << "set yrange [0:]";
    mm_theory1 << "set key on" << "set yrange [0:]";
    mm_experiment
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GMM0"))
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GMM1"))
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
        .Hist(Hist(MC, reaction[0], histpath_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "All")
        .Hist(Hist(MC, reaction[0], histpath_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "Cut1")
        .Hist(Hist(MC, reaction[0], histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "With 2 gammas")
        .Hist(Hist(MC, reaction[0], histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "Cut2")
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
                const hist<> h = Hist(MC, r, histpath_central_reconstr, string("GIM1-Bin-") + to_string(bin_num));
                const hist<> he3 = Hist(MC, r, histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num));
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
                .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM0-Bin-") + to_string(bin_num)), "3He with 2 gammas")
                .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM1-Bin-") + to_string(bin_num)), "Cut2")
                        << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                        << "set xlabel 'gamma-gamma invariant mass, GeV'";
            } else {
                acceptance[i] << point<value<>>(Q, 0.0);
            }
        }
        const hist<> data = Hist(DATA, "C", histpath_central_reconstr, string("GIM1-Bin-") + to_string(bin_num));
        Plot<>(
            Q.Contains(21) ? "He3gg-above-data" : (
                Q.Contains(-39) ? "He3gg-below-data" : (
                    Q.Contains(-3) ? "He3gg-thr-data" : ""
                )
            )
        )
        .Hist(data) << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
        << "set xlabel 'gamma-gamma invariant mass, GeV'";
        ev_am << point<value<>>(Q, value<>::std_error(data.TotalSum().val()));
    }
    Plot<> accplot("He3gg-acceptance");
    accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:]" << "set key on";
    for (size_t i = 0; i < reaction.size(); i++) {
        const auto acc = acceptance[i].YRange(0.0005, INFINITY);
        if (acc.size() > 0) {
            accplot.Hist(acc, reaction[i]);
        }
    }
    const hist<> luminosity = Plotter<>::Instance().GetPoints<value<>>("LUMINOSITYc");
    const hist<> he3etacs = Plotter<>::Instance().GetPoints<value<>>("CS-He3eta-assumed");
    const hist<> known_events =
        luminosity * (runs.first / runs.second)/double(trigger_he3_forward.scaling)
        * (he3etacs * acceptance[1]);
    Plot<>("He3gg-events")
        .Hist(ev_am-known_events, "data-(3He+eta)")
        .Hist(ev_am, "data")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'"<< "set yrange [0:]"
            << "set title '" + runmsg + "'";
}

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
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-6gamma-with-3he");
    vector<string> histpath_reconstr = {"Histograms", "He3nCentralGammas"};
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas6"};
    vector<string> reaction = {"bound1-6g", "He3eta", "He3pi0pi0pi0"};
    hist<> norm = Hist(MC, reaction[0], histpath_reconstr, "0-Reference");
    const auto runs = PresentRuns("C");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";

    Plot gE("He36g-gamma-energy-theory"), gEc("He36g-gamma-energy-theory-cut"),
         gEd("He36g-gamma-energy-data");
    gE << "set key on" << "set title 'Gamma energy. MC'"
       << "set yrange [0:]" << "set xlabel 'gamma energy, GeV'";
    gEc << "set key on" << "set title 'Gamma energy. MC. Cut'"
        << "set yrange [0:]" << "set xlabel 'gamma energy, GeV'";
    gEd << "set key on" << "set title 'Gamma energy. Data'"
        << "set yrange [0:]" << "set xlabel 'gamma energy, GeV'";
    for (size_t i = 0; i < reaction.size(); i++) {
        const auto &r = reaction[i];
        gE.Hist(
            Hist(MC, r, histpath_reconstr, "GammaEnergy")
            / Hist(MC, r, histpath_reconstr, "0-Reference").TotalSum().val()
            , r);
        gEc.Hist(
            Hist(MC, r, histpath_reconstr, "GammaEnergyCut")
            / Hist(MC, r, histpath_reconstr, "0-Reference").TotalSum().val()
            , r);
    }
    gEd.Hist(Hist(DATA, "C", histpath_reconstr, "GammaEnergy")).Hist(Hist(DATA, "C", histpath_reconstr, "GammaEnergyCut"));


    Plot theory("He36g-IMPiDiff0-mc"), theory1("He36g-IMPiDiff1-mc"),
         experiment("He36g-IMPiDiff-data");
    for (const auto &r : reaction) {
        if (r == reaction[0]) {
            Plot("He36g-IMPiDiff-bound-mc")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff3"), "3He 6gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff4"), "3He 3pi^0")
                    << "set key on" << "set yrange [0:]";
        }
        theory.Line(Hist(MC, r, histpath_central_reconstr, "GMMPDiff3").toLine(), r);
        theory1.Line(Hist(MC, r, histpath_central_reconstr, "GMMPDiff4").toLine(), r);
    }
    theory << "set key on" << "set yrange [0:]";
    theory1 << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GMMPDiff3"), "3He 6gamma")
    .Hist(Hist(DATA, "C", histpath_central_reconstr, "GMMPDiff4"), "3He 3pi^0")
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
            Plot(
                Q.Contains(21) ? "He36g-above-he3mm-mc-" + r : (
                    Q.Contains(-39) ? "He36g-below-he3mm-mc-" + r : (
                        Q.Contains(-3) ? "He36g-thr-he3mm-mc-" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(MC, r, histpath_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+6gamma")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He36g-above-ggmm-mc" + r : (
                    Q.Contains(-39) ? "He36g-below-ggmm-mc" + r : (
                        Q.Contains(-3) ? "He36g-thr-ggmm-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM3-Bin-") + to_string(bin_num)), "1 cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "2 cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM5-Bin-") + to_string(bin_num)), "3 cut")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-ggim-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-ggim-mc" + r : (
                        Q.Contains(-3) ? "He3gg-thr-ggim-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM3-Bin-") + to_string(bin_num)), "1 cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num)), "2 cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM5-Bin-") + to_string(bin_num)), "3 cut")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma invariant mass, GeV'";
        }{
            Plot(
                Q.Contains(21) ? "He36g-above-he3mm-data" : (
                    Q.Contains(-39) ? "He36g-below-he3mm-data" : (
                        Q.Contains(-3) ? "He36g-thr-he3mm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(DATA, "C", histpath_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+6gamma")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He36g-above-ggmm-data" : (
                    Q.Contains(-39) ? "He36g-below-ggmm-data" : (
                        Q.Contains(-3) ? "He36g-thr-ggmm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM3-Bin-") + to_string(bin_num)), "1 cut")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "2 cut")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM5-Bin-") + to_string(bin_num)), "3 cut")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He36g-above-ggim-data" : (
                    Q.Contains(-39) ? "He36g-below-ggim-data" : (
                        Q.Contains(-3) ? "He36g-thr-ggim-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM3-Bin-") + to_string(bin_num)), "1 cut")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num)), "2 cut")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM5-Bin-") + to_string(bin_num)), "3 cut")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma invariant mass, GeV'";
        }
        {
            Plot mc_plot(
                Q.Contains(21) ? "He36g-above-im-mc" : (
                    Q.Contains(-39) ? "He36g-below-im-mc" : (
                        Q.Contains(-3) ? "He36g-thr-im-mc" : ""
                    )
                )
            );
            mc_plot << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]"
                    << "set xlabel '6gamma invariant mass, GeV'";
            for (size_t i = 0; i < reaction.size(); i++) {
                const auto &r = reaction[i];
                hist<> Norm = Hist(MC, r, histpath_reconstr, "0-Reference");
                const auto &N = Norm[bin_num].Y();
                if (N.Above(0)) {
                    const hist<> h = Hist(MC, r, histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num));
                    const auto C = std_error(h.TotalSum().val());
                    mc_plot.Hist(h / N, r);
                    acceptance[i] << point<value<>>(Q, C / N);
                } else {
                    acceptance[i] << point<value<>>(Q, 0.0);
                }
            }
        }
        const hist<> data = Hist(DATA, "C", histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num));
        Plot(
            Q.Contains(21) ? "He36g-above-data" : (
                Q.Contains(-39) ? "He36g-below-data" : (
                    Q.Contains(-3) ? "He36g-thr-data" : ""
                )
            )
        )
        .Hist(data)
                << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                << "set xlabel '6gamma invariant mass, GeV'";

        ev_am << point<value<>>(Q, std_error(data.TotalSum().val()));

    }
    Plot accplot("He36g-acceptance");
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
    const hist<> luminosity = Plotter::Instance().GetPoints<value<>>("LUMINOSITYc");
    const hist<> he3etacs = Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed");
    const value<> he3pi0pi0pi0cs={100,25};
    const hist<> he3etaev =
        luminosity * (runs.first / runs.second)
        / double(trigger_he3_forward.scaling)
        * (acceptance[1]*he3etacs);
    const hist<> he3pi0pi0pi0_events =
        luminosity * (runs.first / runs.second)
        / double(trigger_he3_forward.scaling)
        * (acceptance[2]*he3pi0pi0pi0cs);
    Plot("He36g-events")
    .Hist(ev_am, "data")
    .Line(he3etaev.toLine(), "3He+eta estimated")
    .Line(he3pi0pi0pi0_events.toLine(), "3He+3pi^0 estimated")
    .Line(hist<>(he3etaev + he3pi0pi0pi0_events).toLine())
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'";
}

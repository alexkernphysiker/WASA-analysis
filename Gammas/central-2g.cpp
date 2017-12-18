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
    vector<string> histpath_reconstr = {"Histograms", "He3nCentralGammas"};
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound1-2g", "He3eta-gg", "He3pi0pi0", "bound1-6g","He3eta-6g", "He3pi0pi0pi0", "He3pi0"};
    const auto runs = PresentRuns("C");
    const hist<> norm = Hist(MC, reaction[0], histpath_reconstr, "0-Reference");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";

    Plot gE("He3gg-gamma-energy-theory"), gEc("He3gg-gamma-energy-theory-cut"),
         gEd("He3gg-gamma-energy-data");
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
    Plot("He3gg-gamma-count").Hist(Hist(DATA, "C", histpath_reconstr, "GammaCount"))<<"set log y">>"unset log y";

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
                Q.Contains(21) ? "He3gg-above-he3mm-mc-" + r : (
                    Q.Contains(-39) ? "He3gg-below-he3mm-mc-" + r : (
                        Q.Contains(-3) ? "He3gg-thr-he3mm-mc-" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(MC, r, histpath_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM3-Bin-") + to_string(bin_num)), "time cut")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-ggmm-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-ggmm-mc" + r : (
                        Q.Contains(-3) ? "He3gg-thr-ggmm-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM3-Bin-") + to_string(bin_num)), "time cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM5-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-ggim-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-ggim-mc" + r : (
                        Q.Contains(-3) ? "He3gg-thr-ggim-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM3-Bin-") + to_string(bin_num)), "time cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM5-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-tim-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-tim-mc" + r : (
                        Q.Contains(-3) ? "He3gg-thr-tim-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("TIM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("TIM3-Bin-") + to_string(bin_num)), "time cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("TIM5-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'total invariant mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-t-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-t-mc"+r : (
                        Q.Contains(-3) ? "He3gg-thr-t-mc"+r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("t2-Bin-") + to_string(bin_num)))
            .Hist(Hist(MC, r, histpath_central_reconstr, string("t5-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-dt-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-dt-mc"+r : (
                        Q.Contains(-3) ? "He3gg-thr-dt-mc"+r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("dt2-Bin-") + to_string(bin_num)))
            .Hist(Hist(MC, r, histpath_central_reconstr, string("dt5-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";

        }{
            Plot(
                Q.Contains(21) ? "He3gg-above-he3mm-data" : (
                    Q.Contains(-39) ? "He3gg-below-he3mm-data" : (
                        Q.Contains(-3) ? "He3gg-thr-he3mm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(DATA, "C", histpath_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+2gamma")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-ggmm-data" : (
                    Q.Contains(-39) ? "He3gg-below-ggmm-data" : (
                        Q.Contains(-3) ? "He3gg-thr-ggmm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM5-Bin-") + to_string(bin_num)), "MM and IM cuts")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-ggim-data" : (
                    Q.Contains(-39) ? "He3gg-below-ggim-data" : (
                        Q.Contains(-3) ? "He3gg-thr-ggim-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM5-Bin-") + to_string(bin_num)), "MM and IM cuts")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-tim-data" : (
                    Q.Contains(-39) ? "He3gg-below-tim-data" : (
                        Q.Contains(-3) ? "He3gg-thr-tim-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("TIM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("TIM5-Bin-") + to_string(bin_num)), "MM and IM cuts")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'Total invariant mass, GeV'";

            Plot(
                Q.Contains(21) ? "He3gg-above-data-t" : (
                    Q.Contains(-39) ? "He3gg-below-data-t" : (
                        Q.Contains(-3) ? "He3gg-thr-data-t" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("t2-Bin-") + to_string(bin_num)))
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("t5-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-data-dt" : (
                    Q.Contains(-39) ? "He3gg-below-data-dt" : (
                        Q.Contains(-3) ? "He3gg-thr-data-dt" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("dt2-Bin-") + to_string(bin_num)))
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("dt5-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
        }
        {
            Plot mc_plot(
                Q.Contains(21) ? "He3gg-above-im-mc" : (
                    Q.Contains(-39) ? "He3gg-below-im-mc" : (
                        Q.Contains(-3) ? "He3gg-thr-im-mc" : ""
                    )
                )
            );
            Plot mc_plot2(
                Q.Contains(21) ? "He3gg-above-tim-mc" : (
                    Q.Contains(-39) ? "He3gg-below-tim-mc" : (
                        Q.Contains(-3) ? "He3gg-thr-tim-mc" : ""
                    )
                )
            );
            mc_plot << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]"
                    << "set xlabel 'gamma-gamma invariant mass, GeV'";
            mc_plot2 << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'Total invariant mass, GeV'";
            for (size_t i = 0; i < reaction.size(); i++) {
                const auto &r = reaction[i];
                hist<> Norm = Hist(MC, r, histpath_reconstr, "0-Reference");
                const auto &N = Norm[bin_num].Y();
                if (N.Above(0)) {
                    const hist<> h = Hist(MC, r, histpath_central_reconstr, string("GIM5-Bin-") + to_string(bin_num));
                    const auto C = std_error(h.TotalSum().val());
                    mc_plot.Hist(h / N, r);
                    acceptance[i] << point<value<>>(Q, C / N);
                    const hist<> h2 = Hist(MC, r, histpath_central_reconstr, string("TIM5-Bin-") + to_string(bin_num));
                    mc_plot2.Hist(h2 / N, r);
                } else {
                    acceptance[i] << point<value<>>(Q, 0.0);
                }
            }
        }
        const hist<> data = Hist(DATA, "C", histpath_central_reconstr, string("GIM5-Bin-") + to_string(bin_num));
        Plot(
            Q.Contains(21) ? "He3gg-above-data" : (
                Q.Contains(-39) ? "He3gg-below-data" : (
                    Q.Contains(-3) ? "He3gg-thr-data" : ""
                )
            )
        )
        .Hist(data)
                << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                << "set xlabel 'gamma-gamma invariant mass, GeV'";
        ev_am << point<value<>>(Q, std_error(data.TotalSum().val()));
    }
    Plot accplot("He3gg-acceptance");
    accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:]" << "set key on";
    Plot accplot2("He3gg-acceptance2");
    accplot2 << "set title 'Acceptance'"
             << "set xlabel 'Q, MeV'"
             << "set ylabel 'Acceptance, n.d.'"
             << "set yrange [0:0.01]" << "set key on";
    for (size_t i = 0; i < reaction.size(); i++) {
        const auto acc = acceptance[i].YRange(0.000, INFINITY);
        if (acc.size() > 0) {
            accplot.Hist(acc, reaction[i]);
            accplot2.Hist(acc, reaction[i]);
        }
    }

    const hist<> luminosity = Plotter::Instance().GetPoints<value<>>("LUMINOSITYc");
    const hist<> he3etacs = Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed");
    const hist<> known_events =
        luminosity * (runs.first / runs.second)
        / double(trigger_he3_forward.scaling)
        * (acceptance[1]*(he3etacs*0.393));

    Plot("He3gg-events")
    .Hist(ev_am, "data")
    .Hist(known_events, "(3He+eta)_{estimated}")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'";
}

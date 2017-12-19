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

    vector<hist<>> ev_am;
    vector<vector<hist<>>> acceptance;
    const size_t cuts_count=10;
    for(size_t cut_index=0;cut_index<cuts_count;cut_index++){
        acceptance.push_back(vector<hist<>>());
        for (size_t i = 0; i < reaction.size(); i++) {
            acceptance[cut_index].push_back(hist<>());
        }
        ev_am.push_back(hist<>());
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
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "IM and MM cuts")
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
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-tim-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-tim-mc" + r : (
                        Q.Contains(-3) ? "He3gg-thr-tim-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-t-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-t-mc"+r : (
                        Q.Contains(-3) ? "He3gg-thr-t-mc"+r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("t4-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-dt-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-dt-mc"+r : (
                        Q.Contains(-3) ? "He3gg-thr-dt-mc"+r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("dt4-Bin-") + to_string(bin_num)))
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
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "MM and IM cuts")
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
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num)), "MM and IM cuts")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-data" : (
                    Q.Contains(-39) ? "He3gg-below-data" : (
                        Q.Contains(-3) ? "He3gg-thr-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GIM4-Bin-") + to_string(bin_num)))
                << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                << "set xlabel 'gamma-gamma invariant mass, GeV'";
            Plot(
                Q.Contains(21) ? "He3gg-above-tim-data" : (
                    Q.Contains(-39) ? "He3gg-below-tim-data" : (
                        Q.Contains(-3) ? "He3gg-thr-tim-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)), "MM and IM cuts")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";

            Plot(
                Q.Contains(21) ? "He3gg-above-data-t" : (
                    Q.Contains(-39) ? "He3gg-below-data-t" : (
                        Q.Contains(-3) ? "He3gg-thr-data-t" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("t4-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-data-dt" : (
                    Q.Contains(-39) ? "He3gg-below-data-dt" : (
                        Q.Contains(-3) ? "He3gg-thr-data-dt" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("dt4-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
        }
        for(size_t cut_index=0;cut_index<cuts_count;cut_index++){
            const double TIM_DIFF_CUT=0.01*cut_index;
            for (size_t i = 0; i < reaction.size(); i++) {
                const auto &r = reaction[i];
                hist<> Norm = Hist(MC, r, histpath_reconstr, "0-Reference");
                const auto &N = Norm[bin_num].Y();
                if (N.Above(0)) {
                    const hist<> h = Hist(MC, r, histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)).XRange(TIM_DIFF_CUT,INFINITY);
                    acceptance[cut_index][i] << point<value<>>(Q, std_error(h.TotalSum().val()) / N);
                } else {
                    acceptance[cut_index][i] << point<value<>>(Q, 0.0);
                }
            }
            const hist<> data = Hist(DATA, "C", histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)).XRange(TIM_DIFF_CUT,INFINITY);
            ev_am[cut_index] << point<value<>>(Q, std_error(data.TotalSum().val()));
        }
    }
    const hist<> luminosity = Plotter::Instance().GetPoints<value<>>("LUMINOSITYc");
    const hist<> he3etacs = Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed");
    for(size_t cut_index=0;cut_index<cuts_count;cut_index++){
        Plot accplot("He3gg-acceptance-"+to_string(cut_index));
        accplot << "set title 'Cut number "+to_string(cut_index)+"'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0.00001:1]" << "set key on"<<"set log y">>"unset log y";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto acc = acceptance[cut_index][i].YRange(0.00001, INFINITY);
            if (acc.size() > 0) {
                accplot.Hist(acc, reaction[i]);
            }
        }

        const hist<> known_events =
            luminosity * (runs.first / runs.second)
            / double(trigger_he3_forward.scaling)
            * (acceptance[cut_index][1]*(he3etacs*0.393));

        Plot("He3gg-events-"+to_string(cut_index))
        .Hist(ev_am[cut_index], "data")
        .Hist(known_events, "(3He+eta)_{estimated}")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title 'Cut number "+to_string(cut_index)+". "+runmsg+"'";
    }
}

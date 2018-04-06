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
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas6"};
    vector<string> reaction = {"bound3-6g", "He3eta-6g", "He3pi0pi0pi0"};
    hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
    const auto runs = PresentRuns("CC");
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
            Hist(MC, r, histpath_central_reconstr, "GammaEnergy")
            / Hist(MC, r, histpath_central_reconstr, "0-Reference").TotalSum().val()
            , r);
        gEc.Hist(
            Hist(MC, r, histpath_central_reconstr, "GammaEnergyCut")
            / Hist(MC, r, histpath_central_reconstr, "0-Reference").TotalSum().val()
            , r);
    }
    gEd.Hist(Hist(DATA, "CC", histpath_central_reconstr, "GammaEnergy")).Hist(Hist(DATA, "CC", histpath_central_reconstr, "GammaEnergyCut"));


    Plot theory("He36g-IMPiDiff-mc"),experiment("He36g-IMPiDiff-data");
    for (const auto &r : reaction) {
        if (r == reaction[0]) {
            Plot("He36g-IMPiDiff-bound-mc")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff4"))
                    << "set key on" << "set yrange [0:]";
        }
        theory.Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff4"), r);
    }
    theory << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "CC", histpath_central_reconstr, "GMMPDiff4"), "data")
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
            cout<<Qmsg << " plots "<<r<<endl;
            Plot(
                Q.Contains(21) ? "He36g-above-he3mm-mc-" + r : (
                    Q.Contains(-39) ? "He36g-below-he3mm-mc-" + r : (
                        Q.Contains(-3) ? "He36g-thr-he3mm-mc-" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+6gamma")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            cout<<Qmsg << " plots "<<r<<endl;
            Plot(
                Q.Contains(21) ? "He36g-above-ggmm-mc" + r : (
                    Q.Contains(-39) ? "He36g-below-ggmm-mc" + r : (
                        Q.Contains(-3) ? "He36g-thr-ggmm-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "3He+6gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM5-Bin-") + to_string(bin_num)), "6gamma MM cut")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma missing mass, GeV'";
            cout<<Qmsg << " plots "<<r<<endl;
            Plot(
                Q.Contains(21) ? "He36g-above-tim-mc" + r : (
                    Q.Contains(-39) ? "He36g-below-tim-mc" + r : (
                        Q.Contains(-3) ? "He36g-thr-tim-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("TIM5-Bin-") + to_string(bin_num)), "6gamma MM cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("TIM6-Bin-") + to_string(bin_num)), "Total IM cut")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma invariant mass, GeV'";
            cout<<Qmsg << " plots "<<r<<endl;
        }{
            cout<<Qmsg << " plots "<<"data"<<endl;
            Plot(
                Q.Contains(21) ? "He36g-above-he3mm-data" : (
                    Q.Contains(-39) ? "He36g-below-he3mm-data" : (
                        Q.Contains(-3) ? "He36g-thr-he3mm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "CC", histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(DATA, "CC", histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(DATA, "CC", histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+6gamma")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            cout<<Qmsg << " plots "<<"data"<<endl;
            Plot(
                Q.Contains(21) ? "He36g-above-ggmm-data" : (
                    Q.Contains(-39) ? "He36g-below-ggmm-data" : (
                        Q.Contains(-3) ? "He36g-thr-ggmm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "CC", histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "3He+6gamma")
            .Hist(Hist(DATA, "CC", histpath_central_reconstr, string("GMM5-Bin-") + to_string(bin_num)), "6gamma MM cut")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma missing mass, GeV'";
            cout<<Qmsg << " plots "<<"data"<<endl;
        }
        {
            Plot mc_plot(
                Q.Contains(21) ? "He36g-above-tim-mc" : (
                    Q.Contains(-39) ? "He36g-below-tim-mc" : (
                        Q.Contains(-3) ? "He36g-thr-tim-mc" : ""
                    )
                )
            );
            cout<<Qmsg << " plots2 "<<endl;
            mc_plot << "set key on" << "set title '" + Qmsg + ";MC'" << "set yrange [0:]"
                    << "set xlabel '6gamma invariant mass, GeV'";
            for (size_t i = 0; i < reaction.size(); i++) {
                const auto &r = reaction[i];
                cout<<Qmsg << " plots2 "<<r<<endl;
                hist<> Norm = Hist(MC, r, histpath_central_reconstr, "0-Reference");
                const auto &N = Norm[bin_num].Y();
                if (N.Above(0)) {
                    const hist<> h = Hist(MC, r, histpath_central_reconstr, string("TIM6-Bin-") + to_string(bin_num));
                    const auto C = std_error(h.TotalSum().val());
                    mc_plot.Hist(h / N, r);
                    acceptance[i] << point<value<>>(Q, C / N);
                } else {
                    acceptance[i] << point<value<>>(Q, 0.0);
                }
            }
        }
        cout<<Qmsg << " plots3 & events count "<<endl;

        const hist<> DT0 = Hist(DATA, "CC", histpath_central_reconstr, string("dt6-Bin-") + to_string(bin_num)).Scale(50);
        const hist<> DT = Hist(DATA, "CC", histpath_central_reconstr, string("dt7-Bin-") + to_string(bin_num)).Scale(50);
        const hist<> T0 = Hist(DATA, "CC", histpath_central_reconstr, string("t7-Bin-") + to_string(bin_num)).Scale(50);
        const hist<> T = T0.XRange(-10,20);
        const auto TIM=Hist(DATA, "CC", histpath_central_reconstr, string("TIM6-Bin-") + to_string(bin_num));
        cout<<Qmsg << " plots3 & events count "<<endl;

        ev_am<<make_point(Q,std_error(T.TotalSum().val()));
        Plot(
            Q.Contains(21) ? "He36g-above-tim-data" : (
                Q.Contains(-39) ? "He36g-below-tim-data" : (
                    Q.Contains(-3) ? "He36g-thr-tim-data" : ""
                )
            )
        )
            .Hist(TIM)
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma invariant mass, GeV'";
        Plot(
            Q.Contains(21) ? "He36g-above-data1" : (
                Q.Contains(-39) ? "He36g-below-data1" : (
                    Q.Contains(-3) ? "He36g-thr-data1" : ""
                )
            )
        )
        .Hist(DT0).Hist(DT,"Cut")
                << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                << "set xlabel 'gamma-gamma time difference, ns'"<<"set key on";
        Plot(
            Q.Contains(21) ? "He36g-above-data2" : (
                Q.Contains(-39) ? "He36g-below-data2" : (
                    Q.Contains(-3) ? "He36g-thr-data2" : ""
                )
            )
        )
        .Hist(T0).Hist(T,"Cut")
                << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                << "set xlabel '3He-gamma time difference, ns'"<<"set key on";
    }
    cout<<"Final plots"<<endl;
    Plot accplot("He36g-acceptance");
    accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, percents'"
            << "set yrange [0:]" << "set key on";
    for (size_t i = 0; i < reaction.size(); i++) {
        const auto acc = acceptance[i].YRange(0.0001, INFINITY);
        if (acc.size() > 0) {
            accplot.Hist(acc*100, reaction[i]);
        }
    }
    const hist<> luminosity = Plotter::Instance().GetPoints<value<>>("LUMINOSITYc");
    const hist<> luminosity_he = Plotter::Instance().GetPoints<value<>>("LUMINOSITYf");
    const hist<> true_he3eta = luminosity_he
        *hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed"))
            .XRange(luminosity_he.left().X().min(),luminosity_he.right().X().max())
        /trigger_he3_forward.scaling;
    const double branching_ratio=0.32;
    const hist<> known_events = (true_he3eta*branching_ratio)
        *acceptance[1].XRange(true_he3eta.left().X().min(),true_he3eta.right().X().max());
    Plot("He36g-events")
    .Hist(ev_am,"data").Hist(known_events,"3He+eta")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'";
    Plot("He36g-events2")
        .Hist(ev_am.XRange(-50,0),"data, below threshold")
        .Hist(ev_am.XRange(12.5,30)-known_events,"data-3Heeta, upper threshold")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'";
    const auto data_shape=(
        (ev_am*trigger_he3_forward.scaling)/(acceptance[0]*luminosity)
    ).XRange(-45,0);
    Plot("He36g-events-norm-bound")
        .Hist(data_shape,"Data")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]"
            << "set title '"+runmsg+"'";
    Plot("He36g-events-norm2-bound")
        .Hist(data_shape/branching_ratio,"Divided by branching ratio")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]"
            << "set title '"+runmsg+"'";
}

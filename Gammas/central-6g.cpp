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
    const auto runs = PresentRuns("All");
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
    gEd.Hist(Hist(DATA, "All", histpath_central_reconstr, "GammaEnergy")).Hist(Hist(DATA, "All", histpath_central_reconstr, "GammaEnergyCut"));


    Plot theory("He36g-IMPiDiff-mc"),experiment("He36g-IMPiDiff-data");
    for (const auto &r : reaction) {
        if (r == reaction[0]) {
            Plot("He36g-IMPiDiff-bound-mc")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff3"))
                    << "set key on" << "set yrange [0:]";
        }
        theory.Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff3"), r);
    }
    theory << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "All", histpath_central_reconstr, "GMMPDiff4"), "data")
            << "set key on" << "set title '" + runmsg + "'" << "set yrange [0:]";

    hist<> ev_am;
    vector<hist<>> acceptance;
    for (size_t i = 0; i < reaction.size(); i++) {
        acceptance.push_back(hist<>());
    }

    for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            cout<<"All-bins MC plots "<<r<<endl;
            Plot("He36g-he3mm-mc-" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM0"))
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM1"), "cut")
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM3"), "6 gammas required")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass - Q, GeV'";
            Plot("He36g-6gmm-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMM3"), "before cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMM4"), "after cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'";
            Plot("He36g-6gim-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "GIM4"), "before cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GIM5"), "after cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'";
            Plot("He36g-tim-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "TIM5-AllBins"))
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
    }
    {
            Plot("He36g-he3mm-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM0"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM1"), "cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM3"), "6 gammas required")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass - Q, GeV'";
            Plot("He36g-dt-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt2").Scale(2))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            Plot("He36g-t-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t2").Scale(2))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";
            Plot("He36g-6gmm-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GMM3"), "before cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GMM4"), "after cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'";
            Plot("He36g-6gim-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GIM4"), "before cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GIM5"), "after cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass, GeV'";
            Plot("He36g-tim-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "TIM5-AllBins"), "IM and MM cuts")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xrange [-0.3:0.3]" << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
    }
    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
                            << "Q in [" << setprecision(3)
                            << Q.min() << "; " << Q.max() << "] MeV").str();
        cout<<Qmsg << " plots3 & events count "<<endl;
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
                    const hist<> h = Hist(MC, r, histpath_central_reconstr, string("TIM5-Bin-") + to_string(bin_num));
                    const auto C = std_error(h.TotalSum().val());
                    mc_plot.Hist(h / N, r);
                    acceptance[i] << point<value<>>(Q, C / N);
                } else {
                    acceptance[i] << point<value<>>(Q, 0.0);
                }
            }
        }
        cout<<Qmsg << " plots3 & events count "<<endl;
        const auto TIM=Hist(DATA, "All", histpath_central_reconstr, string("TIM5-Bin-") + to_string(bin_num));
        ev_am<<make_point(Q,std_error(TIM.TotalSum().val()));
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

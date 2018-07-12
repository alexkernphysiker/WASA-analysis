// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <math_h/sigma3.h>
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
    vector<string> reaction = {"bound1-6g","bound2-6g","bound3-6g", "He3eta-6g", "He3pi0pi0pi0"};
    hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
    const auto runs = PresentRuns("All");
    const string runmsg = (runs.first==runs.second)?"":"("+to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs)";

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
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff4"))
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff5"),"cut")
                    << "set key on" << "set yrange [0:]";
        }
        theory.Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff5"), r);
    }
    theory << "set key on" << "set yrange [0:]";
    experiment
    .Hist(Hist(DATA, "All", histpath_central_reconstr, "GMMPDiff4"))
    .Hist(Hist(DATA, "All", histpath_central_reconstr, "GMMPDiff5"), "cut")
            << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]";

    hist<> ev_am;
    vector<hist<>> acceptance;
    for (size_t i = 0; i < reaction.size(); i++) {
        acceptance.push_back(hist<>());
    }

    for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            cout<<"All-bins MC plots "<<r<<endl;
            const auto cosb=Hist(MC, r, histpath_central_reconstr,"cosi2")
                        +Hist(MC, r, histpath_central_reconstr,"cosj2")
                        +Hist(MC, r, histpath_central_reconstr,"cosk2");
            const auto cosa=Hist(MC, r, histpath_central_reconstr,"cosi3")
                        +Hist(MC, r, histpath_central_reconstr,"cosj3")
                        +Hist(MC, r, histpath_central_reconstr,"cosk3");
            Plot("He36g-cos-mc"+r)
            .Hist(cosb).Hist(cosa)
                    << "set key on" << "set title '" + r + "'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(gamma-gamma)'";
            Plot("He36g-eta-theta-mc"+r)
            .Hist(Hist(MC, r, histpath_central_reconstr,"ET3"))
            .Hist(Hist(MC, r, histpath_central_reconstr,"ET4"))
                    << "set key on" << "set title '" + r + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'";
            Plot("He36g-he3mm-mc-" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM0"))
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM1"), "cut")
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM2"), "6 gammas required")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass - Q, GeV'"<< "set xrange [0.45:0.57]";
            Plot("He36g-6gmm-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMM5"), "before cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GMM6"), "after cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'"<< "set xrange [2.2:3.4]";
            Plot("He36g-6gim-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "GIM6"), "before cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, "GIM7"), "after cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass + Q, GeV'"<< "set xrange [0.0:1.0]";
            Plot("He36g-tim-mc" + r)
            .Hist(Hist(MC, r, histpath_central_reconstr, "TIM7-AllBins"))
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
            Plot("He36g-he3me-mc"+r)
                .Hist(Hist(MC, r, histpath_central_reconstr, "He3ME7-AllBins"))
                   << "set key on"<<"set yrange [0:]"<< "set title 'Data "+runmsg+"'"
                   << "set xlabel '3He missing energy, GeV'";
    }
    {
            const auto cosb=Hist(DATA, "All", histpath_central_reconstr,"cosi2")
                        +Hist(DATA, "All", histpath_central_reconstr,"cosj2")
                        +Hist(DATA, "All", histpath_central_reconstr,"cosk2");
            const auto cosa=Hist(DATA, "All", histpath_central_reconstr,"cosi3")
                        +Hist(DATA, "All", histpath_central_reconstr,"cosj3")
                        +Hist(DATA, "All", histpath_central_reconstr,"cosk3");
            Plot("He36g-cos-data")
            .Hist(cosb).Hist(cosa)
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(gamma-gamma)'";
            Plot("He36g-eta-theta-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"ET3"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr,"ET4"))
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'";
            Plot("He36g-he3mm-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM0"))
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM1"), "cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM2"), "6 gammas required")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass - Q, GeV'"<< "set xrange [0.45:0.57]";
            Plot("He36g-dt-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt2"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            Plot("He36g-t-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t2"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";
            Plot("He36g-dt-data-final")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt7"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            Plot("He36g-t-data-final")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t7"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";

            Plot("He36g-6gmm-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GMM5"), "before cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GMM6"), "after cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'"<< "set xrange [2.2:3.4]";
            Plot("He36g-6gim-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GIM6"), "before cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "GIM7"), "after cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma invariant mass + Q, GeV'"<< "set xrange [0.0:1.0]";
            Plot("He36g-tim-data")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "TIM7-AllBins"))
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xrange [-0.5:0.5]" << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
            Plot("He36g-he3me-data")
                .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3ME7-AllBins"))
                   << "set key on"<<"set yrange [0:]"<< "set title 'Data "+runmsg+"'"
                   << "set xlabel '3He missing energy, GeV'";
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
                    const hist<> h = Hist(MC, r, histpath_central_reconstr, string("TIM7-Bin-") + to_string(bin_num)).XRange(-0.3,0.3);
                    const auto C = std_error(h.TotalSum().val());
                    mc_plot.Hist(h / N, r);
                    acceptance[i] << point<value<>>(Q, C / N);
                } else {
                    acceptance[i] << point<value<>>(Q, 0.0);
                }
            }
        }
        cout<<Qmsg << " plots3 & events count "<<endl;
        const auto TIM=Hist(DATA, "All", histpath_central_reconstr, string("TIM7-Bin-") + to_string(bin_num)).XRange(-0.3,0.3);
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
    accplot << "set title 'Acceptance'"<<"set key left top"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0:0.1]" << "set key on";
    for (size_t i = 0; i < reaction.size(); i++) {
        const auto acc = acceptance[i].YRange(0.0001, INFINITY);
        if (acc.size() > 0) {
            accplot.Hist(acc, reaction[i]);
        }
    }
    const ext_hist<2> luminosity = Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc");
    const auto luminosity_he = ext_hist<2>(Plotter::Instance().GetPoints<value<>,Uncertainties<2>>("LUMINOSITYc")).XRange(10,30);
    auto true_he3eta = luminosity_he
        *extend_hist<2,2>(hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed")).XRange(10,30))/trigger_he3_forward.scaling;
    while(true_he3eta.left().X().min()>-70.)true_he3eta<<make_point(value<>(true_he3eta.left().X().val()-2.5,2.5),0);

    const auto branching_ratio=uncertainties(0.322,0,0.003);
    const auto known_events = (true_he3eta*branching_ratio)*extend_hist<2,2>(acceptance[3]);
    Plot("He36g-events")
    .Hist(ev_am,"data").Hist_2bars<1,2>(known_events,"3He+eta")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'"<<"set key left top";
    Plot("He36g-events2")
        .Hist_2bars<1,2>(extend_hist<1,2>(ev_am)-known_events,"data-3Heeta")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'";
    const auto data_shape=(extend_hist<1,2>(ev_am)-known_events)*trigger_he3_forward.scaling/(extend_hist<2,2>(acceptance[0])*luminosity);
    Plot("He36g-events-norm-bound")
        .Hist_2bars<1,2>(data_shape.XRange(-70,10),"Data statistical","Data systematic")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]"
            << "set title '"+runmsg+"'"<<"set key right top";
    Plot("He36g-events-norm2-bound")
        .Hist_2bars<1,2>(data_shape.XRange(-70,10)/branching_ratio,"Divided by branching ratio")
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'normalized events amount, nb'" << "set yrange [-10:]"
            << "set title '"+runmsg+"'";
}

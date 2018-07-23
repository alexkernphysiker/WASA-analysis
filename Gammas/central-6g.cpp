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
#include <Parameters/parameters.h>
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

    Plot gE("He36g-gamma-energy-theory",4), gEc("He36g-gamma-energy-theory-cut",4),
         gEd("He36g-gamma-energy-data",4);
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

    {
        Plot theory("He36g-IMPiDiff-mc",4),experiment("He36g-IMPiDiff-data",4);
        for (const auto &r : reaction) {
            const auto imdiff=Hist(MC, r, histpath_central_reconstr, "GMMPDiff4");
            if (r == reaction[0]) {
                Plot("He36g-IMPiDiff-bound-mc",3)
                .Hist(imdiff).Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff5"))
                .Line(Points<>{{getParameter(three_pi0),0.0},{getParameter(three_pi0),hist<>(imdiff.Transponate()).right().X().max()}})
                    << "set key on" << "set yrange [0:]";
            }
            theory.Hist(Hist(MC, r, histpath_central_reconstr, "GMMPDiff5"), r)
            .Line(Points<>{{getParameter(three_pi0),0.0},{getParameter(three_pi0),hist<>(imdiff.Transponate()).right().X().max()}});
        }
        theory << "set key on" << "set yrange [0:]";
        const auto imdiff=Hist(DATA, "All", histpath_central_reconstr, "GMMPDiff4");
        experiment
            .Hist(imdiff).Hist(Hist(DATA, "All", histpath_central_reconstr, "GMMPDiff5"))
            .Line(Points<>{{getParameter(three_pi0),0.0},{getParameter(three_pi0),hist<>(imdiff.Transponate()).right().X().max()}})
                << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]";
    }
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
            Plot("He36g-cos-mc"+r,4)
            .Hist(cosb).Hist(cosa)
                    << "set key on" << "set title '" + r + "'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(gamma-gamma)'";
            const auto he3mm0=Hist(MC, r, histpath_central_reconstr,"He3MM0");
            Plot("He36g-he3mm-mc-" + r,4)
            .Hist(he3mm0).Hist(Hist(MC, r, histpath_central_reconstr,"He3MM1"))
            .Line(Points<>{{getParameter(he3mm_cut),0.0},{getParameter(he3mm_cut),hist<>(he3mm0.Transponate()).right().X().max()}},"cut")
            .Hist(Hist(MC, r, histpath_central_reconstr,"He3MM2"), "6 gammas required")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass - Q, GeV'"<< "set xrange [0.45:0.57]";
            const auto eta_theta=Hist(MC, r, histpath_central_reconstr,"ET3");
            Plot("He36g-eta-theta-mc"+r,4)
            .Hist(eta_theta).Hist(Hist(MC, r, histpath_central_reconstr,"ET4"))
            .Line(Points<>{{getParameter(eta_theta_thr),0.0},{getParameter(eta_theta_thr),hist<>(eta_theta.Transponate()).right().X().max()}},"cut")
                    << "set key on" << "set title '" + r + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'";
            const auto ggmm=Hist(MC, r, histpath_central_reconstr, "GMM5");
            Plot("He36g-6gmm-mc" + r,4)
            .Hist(ggmm).Hist(Hist(MC, r, histpath_central_reconstr, "GMM6"))
            .Line(Points<>{{getParameter(gamma_mm_lo),hist<>(ggmm.Transponate()).right().X().max()},{getParameter(gamma_mm_lo),0.0},
                           {getParameter(gamma_mm_hi),0.0},{getParameter(gamma_mm_hi),hist<>(ggmm.Transponate()).right().X().max()}},"cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '6gamma missing mass, GeV'"<< "set xrange [2.2:3.4]";
            const auto ggim=Hist(MC, r, histpath_central_reconstr, "GIM6");
            Plot("He36g-6gim-mc" + r,4)
            .Hist(ggim).Hist(Hist(MC, r, histpath_central_reconstr, "GIM7"))
            .Line(Points<>{{getParameter(gamma_im_lo6),hist<>(ggim.Transponate()).right().X().max()},{getParameter(gamma_im_lo6),0.0},
                           {getParameter(gamma_im_hi6),0.0},{getParameter(gamma_im_hi6),hist<>(ggim.Transponate()).right().X().max()}},"cut")
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"
                    << "set xlabel '6gamma invariant mass - Q, GeV'"<< "set xrange [0.0:1.0]";
            Plot("He36g-tim-mc" + r,4)
            .Hist(Hist(MC, r, histpath_central_reconstr, "TIM7-AllBins"))
                    << "set key on" << "set title '"+r+"'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'IM(3He+6gamma)-IM(p+d), GeV'";
            Plot("He36g-he3me-mc"+r,4)
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
            Plot("He36g-cos-data",4)
            .Hist(cosb).Hist(cosa)
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-1:1]"
                    << "set xlabel 'cos(gamma-gamma)'";
            const auto he3mm0=Hist(DATA, "All", histpath_central_reconstr, "He3MM0");
            Plot("He36g-he3mm-data",4)
            .Hist(he3mm0).Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM1"))
            .Line(Points<>{{getParameter(he3mm_cut),0.0},{getParameter(he3mm_cut),hist<>(he3mm0.Transponate()).right().X().max()}},"cut")
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "He3MM2"), "6 gammas required")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass - Q, GeV'"<< "set xrange [0.45:0.57]";
            const auto eta_theta=Hist(DATA, "All", histpath_central_reconstr,"ET3");
            Plot("He36g-eta-theta-data",4)
            .Hist(eta_theta).Hist(Hist(DATA, "All", histpath_central_reconstr,"ET4"))
            .Line(Points<>{{getParameter(eta_theta_thr),0.0},{getParameter(eta_theta_thr),hist<>(eta_theta.Transponate()).right().X().max()}},"cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [0:180]"
                    << "set xlabel 'theta(eta) reconstructed'";
            const auto ggmm=Hist(DATA, "All", histpath_central_reconstr, "GMM5");
            Plot("He36g-6gmm-data",4)
            .Hist(ggmm).Hist(Hist(DATA, "All", histpath_central_reconstr, "GMM6"))
            .Line(Points<>{{getParameter(gamma_mm_lo),hist<>(ggmm.Transponate()).right().X().max()},{getParameter(gamma_mm_lo),0.0},
                           {getParameter(gamma_mm_hi),0.0},{getParameter(gamma_mm_hi),hist<>(ggmm.Transponate()).right().X().max()}},"cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma missing mass, GeV'"<< "set xrange [2.2:3.4]";
            const auto ggim=Hist(DATA, "All", histpath_central_reconstr, "GIM6");
            Plot("He36g-6gim-data",4)
            .Hist(ggim).Hist(Hist(DATA, "All", histpath_central_reconstr, "GIM7"))
            .Line(Points<>{{getParameter(gamma_im_lo6),hist<>(ggim.Transponate()).right().X().max()},{getParameter(gamma_im_lo6),0.0},
                           {getParameter(gamma_im_hi6),0.0},{getParameter(gamma_im_hi6),hist<>(ggim.Transponate()).right().X().max()}},"cut")
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma invariant mass - Q, GeV'"<< "set xrange [0.0:1.0]";
            Plot("He36g-tim-data",4)
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "TIM7-AllBins"))
                    << "set key on" << "set title 'Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xrange [-0.5:0.5]" << "set xlabel 'IM(3He+6gamma)-IM(p+d), GeV'";

            const auto DT=Hist(DATA, "All", histpath_central_reconstr, "dt2");
            const auto T=Hist(DATA, "All", histpath_central_reconstr, "t2");
            Plot("He36g-dt-data",4)
            .Hist(DT)
            .Line(Points<>{{getParameter(time_dt),0.0},{getParameter(time_dt),hist<>(DT.Transponate()).right().X().max()}},"cut")
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            Plot("He36g-t-data",4)
            .Hist(T)
            .Line(Points<>{{getParameter(time_t1),hist<>(T.Transponate()).right().X().max()},{getParameter(time_t1),0.0},
                           {getParameter(time_t2),0.0},{getParameter(time_t2),hist<>(T.Transponate()).right().X().max()}},"cut")
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";

            Plot("He36g-dt-data-final",4)
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "dt7"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<< "set key on";
            Plot("He36g-t-data-final",4)
            .Hist(Hist(DATA, "All", histpath_central_reconstr, "t7"))
                    << "set title 'Data " + runmsg + "'"  << "set yrange [0:]"
                    << "set xlabel 'dt 3He-gamma, ns'"<< "set key on";
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
                ),4
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
            ),4
        )
            .Hist(TIM)
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '6gamma invariant mass, GeV'";
    }
    cout<<"Final plots"<<endl;
    Plot accplot("He36g-acceptance",3);
    accplot << "set title 'Acceptance'"<<"set key left top">>"set key right top"
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
    Plot("He36g-events",3)
    .Hist(ev_am,"data").Hist_2bars<1,2>(known_events,"3He+eta")
            <<"set key left top">>"set key right top"
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'"<<"set key left top";
    Plot("He36g-events2",3)
        .Hist_2bars<1,2>(extend_hist<1,2>(ev_am)-known_events,"data-3Heeta")
            <<"set key left top">>"set key right top"
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '" + runmsg + "'";
    const auto data_shape=(extend_hist<1,2>(ev_am)-known_events)*trigger_he3_forward.scaling/(extend_hist<2,2>(acceptance[0])*luminosity);
    Plot("He36g-events-norm-bound",3)
        .Hist_2bars<1,2>(data_shape.XRange(-70,10),"Data statistical","Data systematic")
            <<"set key left top">>"set key right top"
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]"
            << "set title '"+runmsg+"'"<<"set key right top";
    Plot("He36g-events-norm-light",3)
        .Hist(wrap_hist(data_shape).XRange(-70,10))
            <<"set key left top">>"set key right top"
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'normalized events amount, nb'" << "set yrange [0:]"
            << "set title '3He+6gamma"+runmsg+"'"<<"set key right top";
    Plot("He36g-events-norm2-bound",3)
        .Hist_2bars<1,2>(data_shape.XRange(-70,10)/branching_ratio,"Divided by branching ratio")
            <<"set key left top">>"set key right top"
            << "set xlabel 'Q, MeV'" << "set key on"
            << "set ylabel 'normalized events amount, nb'" << "set yrange [-10:]"
            << "set title '"+runmsg+"'";
}

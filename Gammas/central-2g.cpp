// this file is distributed under
// GPL license
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>
#include <gnuplot_wrap.h>
#include <math_h/interpolate.h>
#include <Experiment/experiment_conv.h>
#include <Experiment/str_get.h>
#include <Experiment/gethist.h>
using namespace std;
using namespace ROOT_data;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    Plotter::Instance().SetOutput(ENV(OUTPUT_PLOTS), "central-2gamma");
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound3-2g", "He3eta-gg", "He3pi0pi0", "He3pi0pi0pi0", "He3pi0"};
    const auto runs = PresentRuns("C");
    const hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
    const string runmsg = to_string(int(runs.first)) + " of " + to_string(int(runs.second)) + " runs";
    cout<<"First plots"<<endl;
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
            Hist(MC, r, histpath_central_reconstr, "GammaEnergy")
            / Hist(MC, r, histpath_central_reconstr, "0-Reference").TotalSum().val()
            , r);
        gEc.Hist(
            Hist(MC, r, histpath_central_reconstr, "GammaEnergyCut")
            / Hist(MC, r, histpath_central_reconstr, "0-Reference").TotalSum().val()
            , r);
    }
    gEd.Hist(Hist(DATA, "C", histpath_central_reconstr, "GammaEnergy")).Hist(Hist(DATA, "C", histpath_central_reconstr, "GammaEnergyCut"));
    Plot("He3gg-gamma-count").Hist(Hist(DATA, "C", histpath_central_reconstr, "GammaCount"))<<"set log y">>"unset log y";

    vector<hist<>> ev_am;
    vector<vector<hist<>>> acceptance;
    hist<> he3acc;
    const vector<string> suffix={"-m40","-m20","-0","-20","-40"};
    vector<vector<hist<>>> acc;
    for(size_t i=0;i<suffix.size();i++){
        acc.push_back({});
        for (size_t j = 0; j < reaction.size(); j++)
            acc[i].push_back(hist<>());
        ev_am.push_back(hist<>());
    }

    for (size_t bin_num = 0, bin_count = norm.size(); bin_num < bin_count; bin_num++) {
        const auto &Q = norm[bin_num].X();
        const string Qmsg = static_cast<stringstream &>(stringstream()
                            << "Q in [" << setprecision(3)
                            << Q.min() << "; " << Q.max() << "] MeV").str();
        cout<<Qmsg << " plots"<<endl;
        {
            hist<> he3reg = Hist(MC, reaction[0], histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num));
            he3acc<<make_point(Q,std_error(he3reg.TotalSum().val())/norm[bin_num].Y());
        }
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            cout<<Qmsg << " plots "<<r<<endl;
            Plot(
                Q.Contains(21) ? "He3gg-above-he3mm-mc-" + r : (
                    Q.Contains(-39) ? "He3gg-below-he3mm-mc-" + r : (
                        Q.Contains(-11) ? "He3gg-thr-he3mm-mc-" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+2gamma")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            cout<<Qmsg << " plots "<<r<<endl;
            Plot(
                Q.Contains(21) ? "He3gg-above-ggmm-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-ggmm-mc" + r : (
                        Q.Contains(-11) ? "He3gg-thr-ggmm-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'";
            cout<<Qmsg << " plots "<<r<<endl;
            Plot(
                Q.Contains(21) ? "He3gg-above-tim-mc" + r : (
                    Q.Contains(-39) ? "He3gg-below-tim-mc" + r : (
                        Q.Contains(-11) ? "He3gg-thr-tim-mc" + r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)), "IM and MM cuts")
                    << "set key on" << "set title '" + Qmsg + ";MC " + r + "'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";
            cout<<Qmsg << " plots "<<r<<endl;
            Plot(
                Q.Contains(21) ? "He3gg-above-dt-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-dt-mc"+r : (
                        Q.Contains(-11) ? "He3gg-thr-dt-mc"+r : ""
                    )
                )
            )
            .Hist(Hist(MC, r, histpath_central_reconstr, string("dt4-Bin-") + to_string(bin_num)))
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
            cout<<Qmsg << " plots "<<r<<endl;
        }{
            Plot(
                Q.Contains(21) ? "He3gg-above-he3mm-data" : (
                    Q.Contains(-39) ? "He3gg-below-he3mm-data" : (
                        Q.Contains(-11) ? "He3gg-thr-he3mm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM2-Bin-") + to_string(bin_num)), "3He+2gamma")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '3He missing mass, GeV'";
            cout<<Qmsg << " plots data"<<endl;
            Plot(
                Q.Contains(21) ? "He3gg-above-ggmm-data" : (
                    Q.Contains(-39) ? "He3gg-below-ggmm-data" : (
                        Q.Contains(-11) ? "He3gg-thr-ggmm-data" : ""
                    )
                )
            )
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM2-Bin-") + to_string(bin_num)), "3He+2gamma")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("GMM4-Bin-") + to_string(bin_num)), "MM and IM cuts")
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel '2gamma missing mass, GeV'";
            cout<<Qmsg << " plots data"<<endl;
        }
        for(size_t i = 0; i < reaction.size(); i++) for(size_t a_t=0;a_t<suffix.size();a_t++){
            const auto &r = reaction[i];
            cout<<Qmsg << " plots 2 "<<r<<endl;
            const auto DT=Hist(MC, r, histpath_central_reconstr, "dt5"+to_string(a_t)+"-Bin-"+to_string(bin_num)).Scale(10);
            const auto T=Hist(MC, r, histpath_central_reconstr, "t5"+to_string(a_t)+"-Bin-"+to_string(bin_num)).Scale(10);
            hist<> Norm = Hist(MC, r, histpath_central_reconstr, "0-Reference");
            cout<<Qmsg << " plots 2 "<<r<<endl;
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) acc[a_t][i] << make_point(Q, std_error(DT.TotalSum().val())/N);
            else acc[a_t][i] << make_point(Q, 0.0);
            Plot(
                Q.Contains(21) ? "He3gg-above-dt-final"+suffix[a_t]+"-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-dt-final"+suffix[a_t]+"-mc"+r : (
                        Q.Contains(-11) ? "He3gg-thr-dt-final"+suffix[a_t]+"-mc"+r : ""
                    )
                )
            )
            .Hist(DT)
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-t-final"+suffix[a_t]+"-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-t-final"+suffix[a_t]+"-mc"+r : (
                        Q.Contains(-11) ? "He3gg-thr-t-final"+suffix[a_t]+"-mc"+r : ""
                    )
                )
            )
            .Hist(T)
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
        }
        for(size_t a_t=0;a_t<suffix.size();a_t++){
            cout<<Qmsg<< ";"<<suffix[a_t] << "; plots 3 & events count"<<endl;
            const auto DT0=Hist(DATA, "C", histpath_central_reconstr, "dt5"+to_string(a_t)+"-Bin-"+to_string(bin_num)).Scale(25);
            const auto DT=Hist(DATA, "C", histpath_central_reconstr, "dt6"+to_string(a_t)+"-Bin-"+to_string(bin_num)).Scale(25);
            const auto T0=Hist(DATA, "C", histpath_central_reconstr, "t6"+to_string(a_t)+"-Bin-"+to_string(bin_num)).Scale(25);
            const auto T=T0.XRange(-10,20);
            const auto TIM=Hist(DATA, "C", histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num));
            cout<<Qmsg<< ";"<<suffix[a_t] << "; plots 3 & events count"<<endl;
            ev_am[a_t]<<make_point(Q,std_error(T.TotalSum().val()));
            Plot(
                Q.Contains(21) ? "He3gg-above-tim-data" : (
                    Q.Contains(-39) ? "He3gg-below-tim-data" : (
                        Q.Contains(-11) ? "He3gg-thr-tim-data" : ""
                    )
                )
            )
            .Hist(TIM)
                    << "set key on" << "set title '" + Qmsg + ";Data " + runmsg + "'" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
                    << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";

            Plot(
                Q.Contains(21) ? "He3gg-above-data-dt-final"+suffix[a_t] : (
                    Q.Contains(-39) ? "He3gg-below-data-dt-final"+suffix[a_t] : (
                        Q.Contains(-11) ? "He3gg-thr-data-dt-final"+suffix[a_t] : ""
                    )
                )
            )
            .Hist(DT0).Hist(DT,"Cut")
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'"<<"set key on";
            Plot(
                Q.Contains(21) ? "He3gg-above-data-t-final"+suffix[a_t] : (
                    Q.Contains(-39) ? "He3gg-below-data-t-final"+suffix[a_t] : (
                        Q.Contains(-11) ? "He3gg-thr-data-t-final"+suffix[a_t] : ""
                    )
                )
            )
            .Hist(T0).Hist(T,"Cut")
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'"<<"set key on";
        }

    }
    const hist<> luminosity_he = Plotter::Instance().GetPoints<value<>>("LUMINOSITYf");
    const hist<> true_he3eta = luminosity_he
        *hist<>(Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed"))
            .XRange(luminosity_he.left().X().min(),luminosity_he.right().X().max())
        /trigger_he3_forward.scaling;
    const double branching_ratio=0.39;

    for(size_t a_t=0;a_t<suffix.size();a_t++){
        cout<<suffix[a_t]<< " saving"<<endl;
        const auto TIM=Hist(DATA, "C", histpath_central_reconstr,"TIM4-AllBins");
        const double high=TIM.TransponateAndSort().right().X().max();
        const double cutpos=-0.04+0.02*a_t;
        Plot("He3gg-above-tim-data-allbins"+suffix[a_t]).Hist(TIM)
        .Line({make_point(cutpos,0.),make_point(cutpos,high)})
            << "set key on" << "set yrange [0:]"<< "set xrange [-0.5:0.5]"
            << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d), GeV'";


        Plot accplot("He3gg-acceptance-final"+suffix[a_t]);
        accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, n.d.'"
            << "set yrange [0.0000001:5]" << "set xrange [-70:30]"
            << "set key on"<<"set log y">>"unset log y";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto ac = acc[a_t][i].YRange(0.0000001, INFINITY);
            if (ac.size() > 0) {
                accplot.Hist(ac, reaction[i]);
                Plot("He3gg-acceptance-final"+suffix[a_t]+"-"+reaction[i])
                .Hist(acc[a_t][i],"","He3gg-acceptance"+suffix[a_t]+"-"+to_string(i))
                    << "set title 'Acceptance "+reaction[i]+"'"
                    << "set xlabel 'Q, MeV'"
                    << "set ylabel 'Acceptance, percents'"
                    << "set xrange [-70:30]";
            }
        }
        cout<<suffix[a_t]<< " fitting"<<endl;
        const hist<> known_events = (true_he3eta*branching_ratio)
            *acc[a_t][1].XRange(true_he3eta.left().X().min(),true_he3eta.right().X().max());
        Plot("He3gg-events-final"+suffix[a_t])
            .Hist(ev_am[a_t],"data","He3gg-data"+suffix[a_t])
            .Hist(known_events,"3He+eta")
                << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
                << "set ylabel 'events, n.d.'" << "set yrange [0:]"
                << "set title '"+runmsg+"'";
    }
    cout<<"Final plots"<<endl;
    Plot("He3gg-tube-acc").Hist(he3acc)
            << "set xlabel 'Q, MeV'" << "set xrange [-70:30]"
            << "set ylabel 'Acceptance, n.d.'" << "set yrange [0:1]"
            << "set title 'How many helium ions from mesic nuclei decay would be detected'";
    cout<<"END"<<endl;
}

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
    vector<string> histpath_central_reconstr = {"Histograms", "He3nCentralGammas2"};
    vector<string> reaction = {"bound1-2g", "He3eta-gg", "He3pi0pi0", "He3pi0pi0pi0", "He3pi0"};
    const auto runs = PresentRuns("C");
    const hist<> norm = Hist(MC, reaction[0], histpath_central_reconstr, "0-Reference");
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
    hist<> ev_am1,ev_am2;
    vector<hist<>> acc;
    hist<> e_ev_am1,e_ev_am2;
    vector<hist<>> e_acc;
    const size_t cuts_count=10;
    for(size_t cut_index=0;cut_index<cuts_count;cut_index++){
        acceptance.push_back(vector<hist<>>());
        for (size_t i = 0; i < reaction.size(); i++) {
            acceptance[cut_index].push_back(hist<>());
        }
        ev_am.push_back(hist<>());
    }
    for (size_t i = 0; i < reaction.size(); i++) {
        acc.push_back(hist<>());
        e_acc.push_back(hist<>());
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
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(MC, r, histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
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
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM0-Bin-") + to_string(bin_num)), "3He")
            .Hist(Hist(DATA, "C", histpath_central_reconstr, string("He3MM1-Bin-") + to_string(bin_num)), "3He MM cut")
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
        hist<> ACC,ACCbg,RATIO;
        for(size_t cut_index=0;cut_index<cuts_count;cut_index++){
            const double TIM_DIFF_CUT=0.01*cut_index;
            value<> R=INFINITY;
            for (size_t i = 0; i < reaction.size(); i++) {
                const auto &r = reaction[i];
                hist<> Norm = Hist(MC, r, histpath_central_reconstr, "0-Reference");
                const auto &N = Norm[bin_num].Y();
                if (N.Above(0)) {
                    const hist<> h = Hist(MC, r, histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)).XRange(TIM_DIFF_CUT,INFINITY);
                    const auto a=std_error(h.TotalSum().val()) / N;
                    acceptance[cut_index][i] << point<value<>>(Q,a);
                    if(i==0){
                        R=a;
                        ACC<<make_point(0.001*cut_index,a);
                    }
                    if(i==2){
                        R/=a;
                        ACCbg<<make_point(0.001*cut_index,a);
                    }
                } else {
                    acceptance[cut_index][i] << point<value<>>(Q, 0.0);
                    if(i==0)R=0.0;
                    if(i==2)R=INFINITY;
                }
            }
            if(isfinite(R.val()))RATIO<<make_point(0.001*cut_index,R);
            const hist<> data = Hist(DATA, "C", histpath_central_reconstr, string("TIM4-Bin-") + to_string(bin_num)).XRange(TIM_DIFF_CUT,INFINITY);
            ev_am[cut_index] << point<value<>>(Q, std_error(data.TotalSum().val()));
        }
        if(RATIO.size()>0)
        Plot(
            Q.Contains(21) ? "He3gg-above-acceptance-ratio" : (
                Q.Contains(-39) ? "He3gg-below-acceptance-ratio" : (
                    Q.Contains(-3) ? "He3gg-thr-acceptance-ratio" : ""
                )
            )
        )
        .Hist(RATIO)
                << "set key on" << "set title 'Acceptance ratio" + Qmsg + "'"
                << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, GeV'";
        if(ACC.size()>0)
        Plot(
            Q.Contains(21) ? "He3gg-above-acceptance" : (
                Q.Contains(-39) ? "He3gg-below-acceptance" : (
                    Q.Contains(-3) ? "He3gg-thr-acceptance" : ""
                )
            )
        )
        .Hist(ACC*100)
                << "set key on" << "set title 'Acceptance for foreground " + Qmsg + "'"
                << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, GeV'"
                << "set ylabel 'Acceptance, percents'";
        if(ACCbg.size()>0)
        Plot(
            Q.Contains(21) ? "He3gg-above-acceptance-bg" : (
                Q.Contains(-39) ? "He3gg-below-acceptance-bg" : (
                    Q.Contains(-3) ? "He3gg-thr-acceptance-bg" : ""
                )
            )
        )
        .Hist(ACCbg*100)
                << "set key on" << "set title 'Acceptance for the most intensive background " + Qmsg + "'"
                << "set xlabel 'IM(3He+gamma+gamma)-IM(p+d) cut position, GeV'"
                << "set ylabel 'Acceptance, percents'";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto &r = reaction[i];
            const auto DT=Hist(MC, r, histpath_central_reconstr, string("dt51-Bin-") + to_string(bin_num)).Scale(10);
            const auto T=Hist(MC, r, histpath_central_reconstr, string("t51-Bin-") + to_string(bin_num)).Scale(10);
            const auto e_DT=Hist(MC, r, histpath_central_reconstr, string("dt52-Bin-") + to_string(bin_num)).Scale(10);
            const auto e_T=Hist(MC, r, histpath_central_reconstr, string("t52-Bin-") + to_string(bin_num)).Scale(10);
            hist<> Norm = Hist(MC, r, histpath_central_reconstr, "0-Reference");
            const auto &N = Norm[bin_num].Y();
            if (N.Above(0)) {
                acc[i] << make_point(Q, std_error(DT.TotalSum().val())/N);
                e_acc[i] << make_point(Q, std_error(e_DT.TotalSum().val())/N);
            } else {
                acc[i] << make_point(Q, 0.0);
                e_acc[i] << make_point(Q, 0.0);
            }
            Plot(
                Q.Contains(21) ? "He3gg-above-dt-final-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-dt-final-mc"+r : (
                        Q.Contains(-3) ? "He3gg-thr-dt-final-mc"+r : ""
                    )
                )
            )
            .Hist(DT)
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-t-final-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-t-final-mc"+r : (
                        Q.Contains(-3) ? "He3gg-thr-t-final-mc"+r : ""
                    )
                )
            )
            .Hist(T)
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-dt-final-strict-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-dt-final-strict-mc"+r : (
                        Q.Contains(-3) ? "He3gg-thr-dt-final-strict-mc"+r : ""
                    )
                )
            )
            .Hist(e_DT)
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
            Plot(
                Q.Contains(21) ? "He3gg-above-t-final-strict-mc"+r : (
                    Q.Contains(-39) ? "He3gg-below-t-final-strict-mc"+r : (
                        Q.Contains(-3) ? "He3gg-thr-t-final-strict-mc"+r : ""
                    )
                )
            )
            .Hist(e_T)
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";

        }{
            const auto DT=Hist(DATA, "C", histpath_central_reconstr, string("dt51-Bin-") + to_string(bin_num)).Scale(50);
            const auto T=Hist(DATA, "C", histpath_central_reconstr, string("t51-Bin-") + to_string(bin_num)).Scale(50);
            const auto e_DT=Hist(DATA, "C", histpath_central_reconstr, string("dt52-Bin-") + to_string(bin_num)).Scale(50);
            const auto e_T=Hist(DATA, "C", histpath_central_reconstr, string("t52-Bin-") + to_string(bin_num)).Scale(50);

            const auto DTBG=WeightedAverage<>()<<DT[2].Y()<<DT[3].Y()<<DT[4].Y();
            ev_am1<<make_point(Q,DT[0].Y()+DT[1].Y()-DTBG()*2.0);
            const LinearInterpolation<value<>> TBG=Points<value<>>{T[3],T[6]};
            ev_am2<<make_point(Q,(T[4].Y()-TBG(T[4].X()))+(T[5].Y()-TBG(T[5].X())));
            const hist<> dtbgplot=Points<value<>>{{DT[0].X(),DTBG()},{DT[1].X(),DTBG()}};

            const auto e_DTBG=WeightedAverage<>()<<e_DT[2].Y()<<e_DT[3].Y()<<e_DT[4].Y();
            e_ev_am1<<make_point(Q,e_DT[0].Y()+e_DT[1].Y()-e_DTBG()*2.0);
            const LinearInterpolation<value<>> e_TBG=Points<value<>>{e_T[3],e_T[6]};
            e_ev_am2<<make_point(Q,(e_T[4].Y()-e_TBG(e_T[4].X()))+(e_T[5].Y()-e_TBG(e_T[5].X())));
            const hist<> e_dtbgplot=Points<value<>>{{e_DT[0].X(),e_DTBG()},{e_DT[1].X(),e_DTBG()}};

            Plot(
                Q.Contains(21) ? "He3gg-above-data-dt-final" : (
                    Q.Contains(-39) ? "He3gg-below-data-dt-final" : (
                        Q.Contains(-3) ? "He3gg-thr-data-dt-final" : ""
                    )
                )
            )
            .Hist(DT).Hist(dtbgplot,"background")
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
            const hist<> tbgplot=Points<value<>>{{T[4].X(),TBG(T[4].X())},{T[5].X(),TBG(T[5].X())}};
            Plot(
                Q.Contains(21) ? "He3gg-above-data-t-final" : (
                    Q.Contains(-39) ? "He3gg-below-data-t-final" : (
                        Q.Contains(-3) ? "He3gg-thr-data-t-final" : ""
                    )
                )
            )
            .Hist(T).Hist(tbgplot,"background")
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";

            Plot(
                Q.Contains(21) ? "He3gg-above-data-dt-final-strict" : (
                    Q.Contains(-39) ? "He3gg-below-data-dt-final-strict" : (
                        Q.Contains(-3) ? "He3gg-thr-data-dt-final-strict" : ""
                    )
                )
            )
            .Hist(e_DT).Hist(e_dtbgplot,"background")
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'dt gamma-gamma, ns'";
            const hist<> e_tbgplot=Points<value<>>{{e_T[4].X(),e_TBG(e_T[4].X())},{e_T[5].X(),e_TBG(e_T[5].X())}};
            Plot(
                Q.Contains(21) ? "He3gg-above-data-t-final-strict" : (
                    Q.Contains(-39) ? "He3gg-below-data-t-final-strict" : (
                        Q.Contains(-3) ? "He3gg-thr-data-t-final-strict" : ""
                    )
                )
            )
            .Hist(e_T).Hist(e_tbgplot,"background")
                    << "set title '" + Qmsg + ";" + runmsg + "'" << "set yrange [0:]"
                    << "set xlabel 'quickest gamma - 3he , ns'";
        }

    }
    const hist<> luminosity = Plotter::Instance().GetPoints<value<>>("LUMINOSITYc");
    const hist<> he3etacs = Plotter::Instance().GetPoints<value<>>("CS-He3eta-assumed");
    for(size_t cut_index=0;cut_index<cuts_count;cut_index++){
        Plot accplot("He3gg-acceptance-"+to_string(cut_index));
        accplot << "set title 'Cut number "+to_string(cut_index)+"'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, percents'"
            << "set yrange [0.00001:500]" << "set xrange [-70:30]"
            << "set key on"<<"set log y">>"unset log y";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto acc = acceptance[cut_index][i].YRange(0.0000001, INFINITY);
            if (acc.size() > 0) {
                accplot.Hist(acc*100, reaction[i]);
            }
        }
        Plot("He3gg-events-"+to_string(cut_index))
        .Hist(ev_am[cut_index],"All data")
            << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title 'Cut number "+to_string(cut_index)+". "+runmsg+"'";
    }

    {
        Plot accplot("He3gg-acceptance-final");
        accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, percents'"
            << "set yrange [0.00001:500]" << "set xrange [-70:30]"
            << "set key on"<<"set log y">>"unset log y";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto ac = acc[i].YRange(0.0000001, INFINITY);
            if (ac.size() > 0) {
                accplot.Hist(ac*100, reaction[i]);
                Plot("He3gg-acceptance-final-"+reaction[i]).Hist(ac*100)
                    << "set title 'Acceptance'"
                    << "set xlabel 'Q, MeV'"
                    << "set ylabel 'Acceptance, percents'"
                    << "set xrange [-70:30]";
            }
        }
        const hist<> known_events =
            luminosity
            / double(trigger_he3_forward.scaling)
            * (acc[1]*(he3etacs*0.4));
        Plot("He3gg-events-final").Hist(ev_am1).Hist(ev_am2).Line(known_events.toLine(),"3He+eta")
            << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '"+runmsg+"'";
        Plot("He3gg-events-final-bound").Hist(hist_avr(ev_am1,ev_am2) - known_events)
            << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '"+runmsg+"'";
    }
    {
        Plot accplot("He3gg-acceptance-final-strict");
        accplot << "set title 'Acceptance'"
            << "set xlabel 'Q, MeV'"
            << "set ylabel 'Acceptance, percents'"
            << "set yrange [0.00001:500]" << "set xrange [-70:30]"
            << "set key on"<<"set log y">>"unset log y";
        for (size_t i = 0; i < reaction.size(); i++) {
            const auto ac = e_acc[i].YRange(0.0000001, INFINITY);
            if (ac.size() > 0) {
                accplot.Hist(ac*100, reaction[i]);
                Plot("He3gg-acceptance-final-strict-"+reaction[i]).Hist(ac*100)
                    << "set title 'Acceptance'"
                    << "set xlabel 'Q, MeV'"
                    << "set ylabel 'Acceptance, percents'"
                    << "set xrange [-70:30]";
            }
        }
        const hist<> known_events =
            luminosity
            / double(trigger_he3_forward.scaling)
            * (e_acc[1]*(he3etacs*0.4));
        Plot("He3gg-events-final-strict").Hist(e_ev_am1).Hist(e_ev_am2).Line(known_events.toLine(),"3He+eta")
            << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '"+runmsg+"'";
        Plot("He3gg-events-final-strict-bound").Hist(hist_avr(e_ev_am1,e_ev_am2)-known_events)
            << "set xlabel 'Q, MeV'" << "set key on" << "set xrange [-70:30]"
            << "set ylabel 'events, n.d.'" << "set yrange [0:]"
            << "set title '"+runmsg+"'";
    }
}
